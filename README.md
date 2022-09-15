# Web-APP-test-case-Test-cases-generate-in-parallel
# instruction
code of Web APP test case generations in parallel
# Usage
You need to build the website under test，install the mysql、redis and nginx. you can install phpstudy to Integrate these environments
# Development
python2.7 、 install Levenshtein、selenium、numpy、redis and so on
also need install firefox and geckodriver
# File directory
```
.
##########    the main file
├── Readme.md                        // help
├── dataset                         // Log and results 
├── handle                     
│   ├── async——task.py              // Multi-browser organization scheduling
│   ├── cal_data_fit.py              // Judge the fitness of each path and Gets the coverage of a test case for all paths
│   ├── case.py                      //script of test which can be executed on the browser 
│   ├── generate_seq_ga_check.py               // Legacy algorithm cross-variation and other operations
│   ├── handle_EFSM_module.py             // Process EFSM models
│   ├── main.py                       // Program execution entry
│   ├── recordFun.py                    // record log
│   ├── sensitive_path_info.py           
│   ├── seq_to_script.py                  // script of test which can be executed on the browser 
│   └── obtain_efsm_info2.py              // Gets model information
├── module                               //the txt of model 
└── spath                              //the sensitive path
