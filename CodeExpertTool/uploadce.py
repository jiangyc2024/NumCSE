#dependent: selenium, requests, requests-toolbelt
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import os
from os.path import isfile, isdir, getmtime, basename
import pickle
from time import time, sleep
import requests
import logging
from selenium.common.exceptions import TimeoutException, StaleElementReferenceException
import argparse

def wait_and_find(browser, locator, query_str, timeout = 5):
    return WebDriverWait(browser, timeout, 0.5).until(EC.presence_of_element_located((locator, query_str)))

def wait_and_clickable(browser, locator, query_str, timeout = 5):
    return WebDriverWait(browser, timeout, 0.5).until(EC.element_to_be_clickable((locator, query_str)))

def wait_and_click(browser, locator, query_str, timeout = 5):
    element = WebDriverWait(browser, timeout, 0.5).until(EC.element_to_be_clickable((locator, query_str)))
    browser.execute_script("arguments[0].click();", element)
    return element

def get_project_url(browser, project_name, exer_name):
    browser.get("https://expert.ethz.ch/mycourses/AS20/nmcse/exercises")
    wait_and_click(browser, By.XPATH, f"//tr[./td[2] = '{project_name}']/td[1]/button")
    return wait_and_find(browser, By.XPATH, f"//tr[./td[3] = '{exer_name}']/td[8]//a").get_attribute("href")

def login_codeexpert(browser):
    browser.get("https://expert.ethz.ch/login")
    logging.info("Waiting for user login")
    return WebDriverWait(browser, 9999, 0.5).until(lambda driver: driver.get_cookie("meteor_login_token"))

def login_cached(browser):
    browser.get("https://expert.ethz.ch")
    wait_and_find(browser, By.CSS_SELECTOR, ".ant-btn")
    cache_name = "ce_cookie.pickle"
    if isfile(cache_name):
        mtime = getmtime(cache_name)
        ce_cookie = pickle.load(open(cache_name, "rb"))
        if ce_cookie["expiry"]/1000 + mtime > time():
            browser.add_cookie(ce_cookie)
            return ce_cookie
    ce_cookie = login_codeexpert(browser)
    pickle.dump(ce_cookie, open(cache_name, "wb"))
    return ce_cookie

def upload_file(browser, file_path):
    folder_url = browser.current_url
    exer_key, folder_key = folder_url.split('/')[-1].split('?fileKey=')
    multipart_form_data = (
        ("projectId", (None, exer_key)),
        ("parentKey", (None, folder_key)),
        ("file", (basename(file_path), open(file_path, "rb")))
    )
    login_token = browser.get_cookie("meteor_login_token")
    user_agent = browser.execute_script("return navigator.userAgent;")
    response = requests.post('https://expert.ethz.ch/addToProjectFiles', 
        files=multipart_form_data, 
        cookies={login_token['name']: login_token['value']},
        headers={'Referer': folder_url, 'User-Agent': user_agent})
    logging.info(f"Uploading {file_path}, Response: {response.status_code}")
    return response

def upload_files(browser, files, folder = "."):
    retry = 5
    for _ in range(retry):
        try:
            wait_and_click(browser, By.XPATH, f"//span[@title='{folder}']", timeout=1)
            wait_and_find(browser, By.CSS_SELECTOR, '.ant-upload-btn > input:nth-child(1)', timeout=1)
        except TimeoutException:
            continue
        break
    logging.info(f"Uploading files under {folder}")
    for file in files:
        upload_file(browser, file)
    
def switch_env(browser, is_solution = True, timeout = 20): # True for solution, False for template
    if is_solution:
        wait_and_click(browser, By.XPATH, "//span[text()='Solution']", timeout=timeout)
        wait_and_find(browser, By.XPATH, "//span[text()='Solution' and ../span[contains(@class, 'ant-radio-button-checked')]]", timeout=timeout)
    else:
        wait_and_click(browser, By.XPATH, "//span[text()='Template']", timeout=timeout)
        wait_and_find(browser, By.XPATH, "//span[text()='Template' and ../span[contains(@class, 'ant-radio-button-checked')]]", timeout=timeout)

def delete_file(browser, file_name, retry = 3):
    for _ in range(retry):
        try:
            wait_and_click(browser, By.XPATH, f"//span[@title='{file_name}']", timeout=1)
            wait_and_click(browser, By.XPATH, "//button[@title='delete']", timeout=1)
            wait_and_click(browser, By.CSS_SELECTOR, "button.ant-btn-primary:nth-child(2)", timeout=1)
            while len(browser.find_elements_by_css_selector("button.ant-btn-primary:nth-child(2)")) > 0:
                wait_and_click(browser, By.XPATH, "button.ant-btn-primary:nth-child(2)", timeout=1)
        except TimeoutException:
            continue
        except StaleElementReferenceException:
            refresh_page(browser)
            wait_and_click(browser, By.XPATH, f"//span[@title='{file_name}']", timeout=10)
            continue
        break

def change_file_property(browser, file_name, prop):
    #wait_and_find(browser, By.XPATH, f"//span[@title='{file_name}']")
    #browser.find_elements(By.XPATH, f"//span[@title='{file_name}' and .//svg]")
    prop2dataicon = {"Editable (default file)": "thumbtack",
                     "Editable": "edit",
                     "Read-only": "lock",
                     "Hidden": "eye-slash"}
    wait_and_click(browser, By.XPATH, f"//span[@title='{file_name}']")
    wait_and_find(browser, By.XPATH, f"//span[@title='{file_name}' and contains(@class, 'ant-tree-node-selected')]")
    prop_window = wait_and_clickable(browser, By.CSS_SELECTOR, ".ant-select-selection-item")
    old_property = prop_window.text.strip()
    while prop != old_property:
        prop_window.click()
        wait_and_find(browser, By.CSS_SELECTOR, ".rc-virtual-list-holder-inner")
        #svg not working in XPATH
        wait_and_click(browser, By.XPATH, f"//div[contains(@class, 'ant-select-item-option') and ./div/*[@data-icon = '{prop2dataicon[prop]}']]")
        prop_window = wait_and_clickable(browser, By.CSS_SELECTOR, ".ant-select-selection-item")
        old_property = prop_window.text.strip()

def create_folder(browser, folder_name):
    wait_and_click(browser, By.XPATH, f"//span[@title='.']")
    wait_and_click(browser, By.XPATH, f"//button[@title='add folder']")
    input_n = wait_and_clickable(browser, By.CSS_SELECTOR, "#input")
    input_n.clear()
    input_n.send_keys(folder_name)
    wait_and_click(browser, By.CSS_SELECTOR, "button.ant-btn-primary:nth-child(2)")

def get_files(browser):
    elements = browser.find_elements_by_xpath("//span[contains(@class, 'ant-tree-node-content-wrapper')]")
    return [element.get_attribute("title") for element in elements if element.get_attribute("title") != "."]

def collapse_subfolders(browser):
    open_sign = "//span[contains(@class, 'ant-tree-switcher_open') and ./following-sibling::span[@title != '.']]"
    while browser.find_elements_by_xpath(open_sign):
        try:
            wait_and_click(browser, By.XPATH, open_sign)
        except StaleElementReferenceException:
            pass

def expand_folders(browser):
    close_sign = "//span[contains(@class, 'ant-tree-switcher_close')]"
    while browser.find_elements_by_xpath(close_sign):
        try:
            wait_and_click(browser, By.XPATH, close_sign)
        except StaleElementReferenceException:
            pass

def delete_extra_files(browser, keep_files = []):
    wait_and_click(browser, By.XPATH, "//span[@title='.']")
    expand_folders(browser)
    collapse_subfolders(browser)
    sleep(1)
    files = get_files(browser)
    del_files = set(files) - set(keep_files)
    while del_files:
        for file in del_files:
            delete_file(browser, file)
        files = get_files(browser)
        del_files = set(files) - set(keep_files)

def clear_project(browser):
    delete_extra_files(browser, keep_files=[])

def upload_folder(browser, folder_path, upload_path = ".", file_property = None):
    files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if isfile(os.path.join(folder_path, f))]
    upload_files(browser, files, folder = upload_path)
    for folder in [f for f in os.listdir(folder_path) if isdir(os.path.join(folder_path, f))]:
        create_folder(browser, folder)
        upload_folder(browser, os.path.join(folder_path, folder), upload_path = folder, file_property = file_property)
    sleep(0.5)
    expand_folders(browser)
    if file_property:
        for file in files:
            try:
                change_file_property(browser, basename(file), file_property)
            except TimeoutException:
                continue

def upload_additional_folders(browser, additional_folders, additional_properties):
    files = get_files(browser)
    for folder in additional_folders:
        if (folder not in files) and folder != ".":
            create_folder(browser, folder)
        file_property = None
        if folder in additional_properties:
            file_property = additional_properties[folder]
        upload_folder(browser, additional_folders[folder], upload_path = folder, file_property = file_property)

def get_keep_files(base_folder, additional_folders):
    files = os.listdir(base_folder)
    for key in additional_folders:
        if key == ".":
            files += os.listdir(additional_folders[key])
        else:
            files.append(key)
    return files

def disable_transition(browser):
    browser.execute_script("const styleElement = document.createElement('style');styleElement.setAttribute('id','style-tag');const styleTagCSSes = document.createTextNode('*,:after,:before{-webkit-transition:none!important;-moz-transition:none!important;-ms-transition:none!important;-o-transition:none!important;transition:none!important;-webkit-transform:none!important;-moz-transform:none!important;-ms-transform:none!important;-o-transform:none!important;transform:none!important}');styleElement.appendChild(styleTagCSSes);document.head.appendChild(styleElement);")

def refresh_page(browser):
    browser.refresh()
    disable_transition(browser)

def upload_ceproject(browser, ce_folder_path, additional_folders = {}, additional_properties = {}):
    disable_transition(browser)
    switch_env(browser, False)
    clear_project(browser)
    switch_env(browser, True)
    clear_project(browser)
    solution_path = os.path.join(ce_folder_path, "solution")
    upload_folder(browser, solution_path, upload_path=".")
    upload_additional_folders(browser, additional_folders, additional_properties)
    switch_env(browser, False)
    refresh_page(browser)
    clear_project(browser)
    template_path = os.path.join(ce_folder_path, "template")
    upload_folder(browser, template_path, upload_path=".")
    upload_additional_folders(browser, additional_folders, additional_properties)
    refresh_page(browser)
    switch_env(browser, True)
    delete_extra_files(browser, keep_files=get_keep_files(solution_path, additional_folders))
    switch_env(browser, False)
    delete_extra_files(browser, keep_files=get_keep_files(template_path, additional_folders))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Code Expert Upload Tool')
    parser.add_argument("-d", "--driver_path", help="Executable path of geckodriver")
    parser.add_argument("-u", "--exercise_url", help="Url of target Code Expert exercise")
    parser.add_argument("-p", "--exercise_path", help="Path of exercise code (the path containing `solution` and `template`)")
    parser.add_argument("-t", "--testscript_path", help="Path of testing scripts (normally is Testing folder under NumCSE repository")
    args = parser.parse_args()
    additional_folders = {".": args.testscript_path}
    additional_properties = {".": "Hidden"}
    browser = webdriver.Firefox(executable_path = args.driver_path)
    login_codeexpert(browser)
    browser.get(args.exercise_url)
    upload_ceproject(browser, args.exercise_path, additional_folders, additional_properties)