# Code Expert upload tool

This is a python script that helps uploading CE projects.

## Dependency
Please install following dependent packages in python (>=3.6) using `pip install --user` or `conda install `
- selenium
- requests

Selenium also requires a Firefox browser driver. You need to install the Firefox browser 
and download the [geckodriver](https://github.com/mozilla/geckodriver/releases). Please 
note down the executable path of `geckodriver` and specify it to this script.

## Get Started
1. Login Code Expert and create an empty exercise
2. Click into the exercise and copy its url, should be in the form of "https://expert.ethz.ch/ide2/<HASHCODE>"
3. Find the absolute path for the exercise folder, should include `solution` and `template` subfolders
4. Find the absolute path for Testing script folder, should be `Testing` folder inside `NumCSE` repo
5. Call this script with infomation found above: `python uploadce.py --driver_path <DRIVER_PATH> --exercise_url <EXERCISE_URL> --exercise_path <EXERCISE_PATH> --testscript_path <TESTSCRIPT_PATH>`
6. The script will open an automated Firefox window, and direct to Code Expert login page, please login with your ETH identity and wait for the uploading process

## Problems
- When having timeout errors, rerunning the script may solve the issue
- This script cannot setup the visibility properties for all the files, you may need to config manually for some files
