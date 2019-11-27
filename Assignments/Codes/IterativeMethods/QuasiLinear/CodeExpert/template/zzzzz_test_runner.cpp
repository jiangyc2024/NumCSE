/*
The MIT License (MIT)

Copyright (c) 2016-2018 ETH Zurich

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
Current Version 1.12

Changes since 1.12:
- remove HTML output
Changes since 1.11:
- close file descriptor to feed input data after all data has been sent
Changes since 1.10:
- output score only in percentage similar to code expert
Changes since 1.9:
- added test for submission -- submissions automagically run the tests
- Changed for C++11 and pedantic warnings compliance
Changes since 1.8:
- Randomize the first 64KB of stack before test execution to trigger bugs due
  to uninitialized variables.
Changes since 1.7:
- Print score in test results
Changes since 1.6:
- Append all characters the program under test outputs, i.e., do not stop if a
  zero character is found.
- Replace non-printable characters (<#32) except the line feed (#10), with a
  reversed question mark.
Changes since 1.5:
- Set score to 1.0 if total score is zero, e.g., if no test cases are defined.
Changes since 1.4:
- Switched from non blocking I/O to polling I/O. This to avoid a race condition
  between program under test and test runner. With non-blocking I/O the program
  under test may read faster from stdin than the test runner can provide input.
  Thus it can miss certain inputs and provide random results.
- Clear output before program execution.
- Fixed some warnings
Changes since 1.3:
- Show a disclaimer after test results.
Changes since 1.2:
- Added automatic floating point value detection. If both, the expected and the
  actual output contain a floating point value it is compared using absolute
  and relative precision (default epsilon - 1e-6). In the expected output,
  floating point values are distinguished from integer values by the means of
  the decimal dot or the exponent.
Changes since 1.1:
- Non testing mode: Print result string after for NO_TEST_TIMEOUT that lies just
  under the codeboard timeout of 20s without aborting the program as the user
  might be providing input.

Usage:
- Tests are defined in the file tests.csv. They are defined using an array of
  C-strings that is included upon compilation. A test case consists of (in this order):
  - Name
  - Input
  - Expected output
  - Score
  - Timeout
*/

// standard C/C++ includes
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cerrno>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// linux specific includes
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <poll.h>

// usec
#define SLEEP_MIN 1
#define SLEEP_MAX 50
#define SLEEP_INC 1

#define NO_TEST_TIMEOUT 15000

// required floating point precision (absolute and relative).
#define EPS 1e-9


int main();

namespace test_runner {


enum Action {
  RUN = 0,
  TEST = 1,
  SUBMIT = 2
};

Action GetAction(){
  std::string action = std::getenv("ACTION");
  if (action == "test") return TEST;
  else if (action == "submit") return SUBMIT;
  else return RUN;
}

enum TestStatus {
    STATUS_OK       = 0,
    STATUS_WRONG    = 1,
    STATUS_ERR      = 2,
    STATUS_TIMEOUT  = 3,
    STATUS_INTERNAL = 4,    
};


const auto red = "\033[31;1m";
const auto indent = "    ";
const auto green = "\033[32;1m";
const auto yellow = "\033[33;1m";
const auto reset = "\033[0m";
const auto sep = "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━";

// caution: next line is difficult to edit
const std::vector<std::string> blocks = {"", "▏","▎","▍","▌","▋","▊","▉"};  

struct TestCase {
    // static data
    std::string name;
    std::string in;
    std::string expect_out;
    int timeout;
    int score;
    
    // runtime data
    std::string actual_out;
    TestStatus result;
    std::string result_msg;
    int running_time;
    
    TestCase(std::string _name, std::string _in, std::string _expect_out, int _score, int _timeout) :
        name(_name), in(_in), expect_out(_expect_out), timeout(_timeout), score(_score) {}
        
    // Sets internal error status, and adds a detail message to the actual
    // output.
    void set_status(const char* message, int retval, int err);
    
    // Returns true if status indicates success.
    bool has_passed();

    // Returns the result of a test as string.
    std::string get_result_string();

    // Returns the time since the epoch in milliseconds.
    static long long int current_millis();

    // Checks whether the result of the testcase is valid. The comparison is
    // done on a per line basis. Superfluous content and whitespace, empty lines,
    // and, casing are ignored.
    static bool is_valid_result(const std::string& expected, const std::string& actual);
    
    // Return true, if all elements from expected_line are contained in
    // actual_line. Elements are seprated by whitespace. Superfluous whitespace
    // is ignored.
    static bool contains_elements(const std::string& expected_line, const std::string& actual_line);

    // Return true if a string contains only whitespace.
    static bool is_empty(const std::string& s);
    
    // Returns true if two strings are considered equal.
    // This function iqnores the case of alphabetic characters.
    static bool is_equal(const std::string& s1, const std::string& s2);

    // Returns true, if the string could be converted to a double. No leading or
    // trailing whitespace or characters must be present as otherwise conversion
    // fails.
    static bool to_double(double& d, const std::string& s);
    
    // Returns true, if two floating point values are considered equal. They are
    // considered equal if either the absolute or relative difference lies below
    // a threshold (eps = 1e-6).
    static bool is_double_equal(double d1, double d2);

    // Monitors the execution of the program under test.
    void monitor_child_execution(pid_t pid, int write_fd, int read_fd);
    
    // Update actual output with contents of buffer 'buf' having length 'len'.
    void upd_actual_out(const char buf[], const int len);
    
    // Executes the test case.
    void execute();
};


class TestRunner {
    std::vector<TestCase> test_cases;

    // "Loads" the test cases. As code board currently does not support
    // loading from existing file we include the data in the test runner.
    void load_test_cases();

    // Executes all test cases.
    void execute_test_cases();

    // Prints the result of a test cases in HTML format.
    void print_results();

    // Exits the program immediately.
    static void shutdown(int rc);

    // Replace non-printable characters in the string with HTML codes.
    std::string replace_non_printable(const std::string &s);
    
    // fill the first 64KB of the stack with random numbers
    static void randomize_stack();
public:
    // Starts the benchmark runner.
    void run();

    // Executes the main function of the porgram under test, and exits the
    // process.
    static void execute_child();
    
    TestRunner();
};

void TestCase::set_status(const char* message, int retval, int err) {
    std::stringstream out;
    out << message << strerror(err) << "(" << err << "), return value: " << retval << "\n";
    result_msg = out.str();
    result = STATUS_INTERNAL;
}

bool TestCase::has_passed() {
    return result == STATUS_OK;
}

std::string TestCase::get_result_string() {
    switch(result) {
        case STATUS_OK:       return std::string("passed");        break;
        case STATUS_WRONG:    return std::string("wrong answer");   break;
        case STATUS_ERR:      return std::string("program error");  break;
        case STATUS_TIMEOUT:  return std::string("timeout");        break;
        case STATUS_INTERNAL: return std::string("internal error"); break;
        default:              return std::string("unknown");      break;
    }
}

long long int TestCase::current_millis() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long int) tv.tv_sec * 1000 + (long long int) tv.tv_usec / 1000;
}

bool TestCase::is_valid_result(const std::string& expected, const std::string& actual) {
    std::stringstream expected_stream(expected);
    std::stringstream actual_stream(actual);

    bool result = true;
    do {
        // get next line while skipping empty lines
        std::string expected_line;
        do {
            getline(expected_stream, expected_line);
        } while (!expected_stream.eof() && is_empty(expected_line));
        
        //std::cout << "search line: " << expected_line << "\n";
        
        // get next actual line matching the expected line or stop
        if (!is_empty(expected_line)) {
            std::string actual_line;
            do {
                getline(actual_stream, actual_line);
            } while (!actual_stream.eof() && !contains_elements(expected_line, actual_line));
            
            result &= contains_elements(expected_line, actual_line);
            
            //std::cout << "found: " << actual_line << " / status: " << result << "\n";
        }
        
    } while(!expected_stream.eof());
    
    return result;
}

bool TestCase::is_empty(const std::string& s) {
    return s.find_first_not_of("\t\n ") == std::string::npos;
}
    
// checks whether all elements from expected_line are contained in actual_line
// whitespace is ignored
bool TestCase::contains_elements(const std::string& expected_line, const std::string& actual_line) {
    std::stringstream expect_stream(expected_line);
    std::stringstream actual_stream(actual_line);
    std::string expect;
    std::string actual;
    expect_stream >> expect;
    actual_stream >> actual;
    
    // find first element
    while (!is_equal(expect, actual) && expect_stream && actual_stream) {
        actual_stream >> actual;
    }
    
    // if found then find remaining elements
    bool result;
    if (is_equal(expect, actual)) {
        while(is_equal(expect, actual) && expect_stream && actual_stream) {
            expect_stream >> expect;
            actual_stream >> actual;
        }
        result = is_equal(expect, actual) || (!expect_stream && actual_stream);
    } else {
        result = false;
    }
    
    return result;
}

bool TestCase::to_double(double& d, const std::string& s) {
    const char* c_str = s.c_str();
    char* next;
    d = strtod(c_str, &next); // next: position after float
    bool success = (unsigned)(next - c_str) == s.length(); // do not accept partial floats
    return success;
}

bool TestCase::is_double_equal(double d1, double d2) {
    bool is_equal;
    
    if (std::isnan(d1) || std::isnan(d2)) {
        // if x is NaN then x == x doesn't hold.
        is_equal = std::isnan(d1) && std::isnan(d2);
    } else if (std::isinf(d1) || std::isinf(d2)) {
        is_equal = d1 == d2;
    } else {
        is_equal = fabs(d1 - d2) <= EPS || fabs(d1 - d2) / fabs(d1) <= EPS;
    }
    
    return is_equal;
}

bool TestCase::is_equal(const std::string& s1, const std::string& s2) {
    std::string lower_s1(s1);
    std::string lower_s2(s2);
    std::transform(lower_s1.begin(), lower_s1.end(), lower_s1.begin(), ::tolower);
    std::transform(lower_s2.begin(), lower_s2.end(), lower_s2.begin(), ::tolower);
    
    double d1, d2;
    bool is_fp = s1.find_first_of(".eE") != std::string::npos;
    if (is_fp && to_double(d1, s1) && to_double(d2, s2)) {
        return is_double_equal(d1, d2);
    } else {
        return lower_s1 == lower_s2;
    }
}

void TestCase::monitor_child_execution(pid_t pid, int write_fd, int read_fd) {
    long long int start_time = current_millis();

    const char* in_cstr = in.c_str();
    int in_pos          = 0;
    int in_remn         = in.size();

    int poll_timeout = SLEEP_MIN;
    
    pollfd watched_fds[2];
    watched_fds[0].fd     = write_fd;
    watched_fds[0].events = POLLOUT;
    watched_fds[1].fd     = read_fd;
    watched_fds[1].events = POLLIN;

    // disable buffering of input
    FILE* fd = fdopen(write_fd, "w");
    setvbuf(fd, NULL, _IONBF, 0);

    bool child_running = true;
    
    while(true) {
        int poll_ret = poll(watched_fds, 2, poll_timeout);
        if (poll_ret == -1) {
            set_status("poll operation failed: ", poll_ret, errno);
            break;
        }
        
        // write test data to child
        if ((watched_fds[0].revents & POLLOUT) && in_remn > 0) {
            int res = write(write_fd, in_cstr + in_pos, in_remn);
            if (res == -1) {
                if (errno != EBADF) {
                    int err = errno;
                    set_status("write operation failed: ", res, err);
                    break;
                }
            } else {
                in_pos  += res;
                in_remn -= res;
            }
        }
        
        // close write_fd if not yet closed and no more data to write
        if (in_remn == 0) {
            close(write_fd);
            in_remn = -1;
        }

        // read output of child
        char buffer[4096];
        int count = 0;
        if ((watched_fds[1].revents & POLLIN)) {
            count = read(read_fd, buffer, sizeof(buffer) - 1);
            if (count > 0) {
                upd_actual_out(buffer, count);

            } else if (count < 0) {
                int err = errno;
                if (err != EBADF) {
                    set_status("read operation failed: ", count, err);
                    break;
                }
            }
        }
        
        // exit after all output has been read
        if (count == 0 && !child_running) {
            break;
        }
        
        // evaluate child status
        if (child_running) {
            int wstatus = 0;
            int child_pid = waitpid(pid, &wstatus, WNOHANG);
            running_time = current_millis() - start_time;
            if (child_pid > 0) {
                if (WIFEXITED(wstatus)) {
                    int rc = WEXITSTATUS(wstatus);
                    if (rc == 0) {
                        result = STATUS_OK;
                    } else {
                        result_msg = "error: program terminated with non-zero exit code";
                        result = STATUS_ERR;
                    
                    }
                } else if (WIFSIGNALED(wstatus)) {
                    int signal = WTERMSIG(wstatus);
                    result_msg  = "program terminated by signal: ";
                    result_msg += strsignal(signal);
                    result = STATUS_ERR;
                    
                } else {
                    result_msg = "error: program terminated abnormally (unknown reason)";
                    result = STATUS_ERR;
                    
                }
                child_running = false;
            
            } else if (child_pid == 0) {
                // check for timeout
                if (running_time >= timeout) {
                    kill(pid, SIGKILL);
                    result = STATUS_TIMEOUT;
                    child_running = false;
                }
            } else {
                set_status("error while waiting for child process: ", child_pid, errno);
                break;
            }
        }
        
        // increase sleep time with each execution until SLEEP_MAX is reached.
        // reduces load of monitor function on test system, while retaining
        // a more or less accurate time measurement.
        if (poll_timeout < SLEEP_MAX) {
            poll_timeout += SLEEP_INC;
        }
    }

    // compute result
    if (result == STATUS_OK) {
        if (!is_valid_result(expect_out, actual_out)) {
            result = STATUS_WRONG;
        }
    }
    if (result_msg != "") {
        actual_out += result_msg;
    }
}

void TestCase::upd_actual_out(const char buf[], const int len) {
    for (int i = 0; i < len; ++i) {
        actual_out += buf[i];
    }
}

void TestCase::execute() {
    #define NUM_PIPES          2
    #define PARENT_WRITE_PIPE  0
    #define PARENT_READ_PIPE   1
    
    // pipes for parent to write and read
    int pipes[NUM_PIPES][2];
    pipe(pipes[PARENT_READ_PIPE]);
    pipe(pipes[PARENT_WRITE_PIPE]);

    #define READ_FD  0
    #define WRITE_FD 1

    #define PARENT_READ_FD  ( pipes[PARENT_READ_PIPE][READ_FD]   )
    #define PARENT_WRITE_FD ( pipes[PARENT_WRITE_PIPE][WRITE_FD] )

    #define CHILD_READ_FD   ( pipes[PARENT_WRITE_PIPE][READ_FD]  )
    #define CHILD_WRITE_FD  ( pipes[PARENT_READ_PIPE][WRITE_FD]  )

    // flush stdout before forking,
    // otherwise child inherits buffered output of parent
    fflush(stdout);
    pid_t pid = fork();
    
    if (pid > 0) {
        // PARENT PROCESS
        
        // close fds not required by parent
        close(CHILD_READ_FD);
        close(CHILD_WRITE_FD);
        
        monitor_child_execution(pid, PARENT_WRITE_FD, PARENT_READ_FD);
        
    } else if (pid == 0) {
        // CHILD PROCESS
        
        // replace stdin and stdout
        dup2(CHILD_READ_FD, STDIN_FILENO);
        dup2(CHILD_WRITE_FD, STDOUT_FILENO);

        // close fds not required by child
        // we don't want the executed program to know these existed
        close(CHILD_READ_FD);
        close(CHILD_WRITE_FD);
        close(PARENT_READ_FD);
        close(PARENT_WRITE_FD);

        TestRunner::execute_child(); // exits

    } else {
        // close fds not required by parent
        close(CHILD_READ_FD);
        close(CHILD_WRITE_FD);
        
        set_status("failed to fork child process: ", pid, errno);
    }
    
    close(PARENT_READ_FD);
    close(PARENT_WRITE_FD);
}


void TestRunner::load_test_cases() {
    const char *raw_data[] = {
       #include "tests.csv"
       NULL
    };
    
    int idx = 0;
    while(raw_data[idx] != NULL) {
        std::string name       = raw_data[idx++];
        std::string in         = raw_data[idx++];
        in += "\n";
        std::string expect_out = raw_data[idx++];
        int score              = strtol(raw_data[idx++], NULL, 10);
        int timeout            = strtol(raw_data[idx++], NULL, 10);
        
        TestCase test_case(name, in, expect_out, score, timeout);
        test_cases.push_back(test_case);
    }
}

void TestRunner::execute_test_cases() {
    std::cout << "Running tests";

    for(std::vector<TestCase>::iterator it = test_cases.begin(); it != test_cases.end(); ++it) {
        (*it).execute();
        std::cout << '.';
    }
    std::cout << "\n\n";
}

std::string newLn(const std::string& str)
{
    if (str.empty() || str.back() != '\n')
      return "\n";
    else
      return "";
}

std::string progressBar(double val){
  std::string result ="[";
  const int characters = 20; 
  const int maxticks = characters * 8;
  int i=0;
  int ticks = round(val * maxticks); // 8 characters
  result += green;
  for (; i< ticks; i+=8){
    result += blocks[7];
  }
  result += blocks[ticks % 8];
  result += reset;
  for (;i<maxticks; i+=8){
    result += " ";
  }
  result += "]";
  return result;
}

void TestRunner::print_results() {
    
    int total           = 0;
    int succeeded       = 0;
    int total_score     = 0;
    int succeeded_score = 0;
    
    for(std::vector<TestCase>::iterator it = test_cases.begin(); it != test_cases.end(); ++it) {
        TestCase& test_case = *it;
        if (test_case.result == STATUS_OK) {
            ++succeeded;
            succeeded_score += test_case.score;
        }
        ++total;
        total_score += test_case.score;
    }

    double score;
    if (total_score != 0) {
        score = (double) succeeded_score / total_score;
    } else {
        score = 1.0;
    }
    
    // print test cases
    for(std::vector<TestCase>::iterator it = test_cases.begin(); it != test_cases.end(); ++it) {
        TestCase& test_case = *it;
        
        std::string result   = test_case.get_result_string();
        
        if (test_case.has_passed()){
          std::cout << test_case.name << ": ";
          std::cout << green << result << reset << std::endl;
        }
        
    }
    
    std::cout << sep << std::endl;
    
    // print test cases
    for(std::vector<TestCase>::iterator it = test_cases.begin(); it != test_cases.end(); ++it) {
        TestCase& test_case = *it;
        
        std::string result   = test_case.get_result_string();
        
        if (!test_case.has_passed()){
          std::cout << test_case.name << ": ";
          std::cout << red << result << reset << std::endl;
          std::cout << indent << yellow << "input:" << reset << std::endl;
          std::cout << indent << test_case.in << newLn(test_case.in);
          std::cout << indent << yellow << "expected output:" << reset << std::endl;
          std::cout << indent << test_case.expect_out <<  newLn(test_case.expect_out);
          std::cout << indent << yellow << "actual output:" << reset << std::endl;
          std::cout << indent << test_case.actual_out <<  newLn(test_case.actual_out);
          std::cout << sep << std::endl;
        }
        
    }
    
  std::cout <<  "\nTests result: passed " << succeeded << " of " << total << " / score: " << round(score * 100.0) << "% " << progressBar(score) << std::endl;

    
  if (GetAction()==SUBMIT){
    std::cout << std::setprecision(4) << "<cx:result value=\"" << score << "\" />";
  }
}

std::string TestRunner::replace_non_printable(const std::string &s) {
	std::string out;
    for (unsigned int i = 0; i < s.length(); ++i) {
        char c = s[i];
        if (c > 31 /*print*/ || c == 10 /*LF*/) {
            out += c; // append to output            
        } else {
            out += "&#191;"; // reversed question mark
        }
    }
	return out;
}


void TestRunner::randomize_stack() {
    std::srand(std::time(0)); // use current time as seed for random generator
    int mem[65536];
    for(int i = 0; i < 65536; ++i) {
        mem[i] = std::rand();
    }
    if (mem[10] == 0) mem[10] = 1; // to avoid triggering the unused variable warning
}

void TestRunner::execute_child() {
    // do not buffer standard out, this provides output, even
    // if the child process is killed.
    setvbuf(stdout, NULL, _IONBF, 0);
    
    // randomize the first 64KB of the stack to enable repeatable 
    // errors in case of uninitialized variables
    randomize_stack();
    
    // __extension__ used in order to disable a warning ISO C++ forbids taking address of function '::main'
    __extension__ int rc = main();
    
    exit(rc);
}

void TestRunner::shutdown(int rc) {
    // flush manually as _exit does not necessarily flush according to specification
    fflush(stdout);
    fflush(stderr);

    // Use _exit to disallow overwriting of result string.
    // Exits without calling any:
    // - at_exit handlers
    // - at_quick_exit handlers
    // - static destructors
    _exit(rc);
}

void TestRunner::run() {
    load_test_cases();
    execute_test_cases();
    print_results();
    shutdown(0);
}

TestRunner::TestRunner() {
    int rc;

    //std::cout << TestCase::contains_elements("Hello World: Heinz ", "hello world: Hein") << "\n";
    //std::cout << TestCase::is_valid_result("\nHello World: Heinz\nz\n", "\nhello world: Heinz\nabc z\n") << "\n";
    
    Action action = GetAction();
    if (action == TEST || action == SUBMIT) {
        run();
        rc = 0;
    } else {
        // execute directly
        pid_t pid = fork();
        if (pid > 0) {
            // parent process
            long long int start_time = TestCase::current_millis();
            long long int running_time = 0;
            
            // ensure that the is a result string returned by printing it
            // after a predetermined timeout
            while(true) {
                running_time = TestCase::current_millis() - start_time;
                
                // check if child is still running
                int wstatus = 0;
                pid_t child_pid;
                
                if (running_time < NO_TEST_TIMEOUT) {
                    child_pid = waitpid(pid, &wstatus, WNOHANG);
                } else {
                    // before timeout occurs print result string, then wait indefinitly
                    if (GetAction() == SUBMIT){
                      std:: cout << "<cx:result value=\"0\" />" << std::endl;
                    }
                    fflush(stdout);
                    child_pid = waitpid(pid, &wstatus, 0);
                }

                if (child_pid > 0) {
                    if (WIFEXITED(wstatus)) {
                        rc = WEXITSTATUS(wstatus);
                        
                    } else if (WIFSIGNALED(wstatus)) {
                        int signal = WTERMSIG(wstatus);
                        std::cerr << "program terminated by signal: " << strsignal(signal) << "\n";
                        rc = 128 + signal;
                        
                    } else {
                        std::cerr << "error: program terminated abnormally (unknown reason)\n";
                        rc = -1;
                    }
                    break;
            
                } else if (child_pid < 0) {
                    std::cerr << "error: waiting for child process to finish failed\n";
                    rc = -1;
                    break;
                }
                
                usleep(SLEEP_MAX * 1000); // sleep to reduce load on test system
            }
            if (GetAction() == SUBMIT){
              std:: cout << "<cx:result value=\"0\" />" << std::endl;
            }

        } else if (pid == 0) {
            test_runner::TestRunner::execute_child();
            
        } else {
            // failed to start child
            std::cerr << "error: failed to fork child process\n";
            rc = -1;
        }
    }
    
    TestRunner::shutdown(rc);
}

// execute the test runner just before main starts
TestRunner run;

} // namespace test_runner
