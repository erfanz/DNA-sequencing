# DNA-sequencing

About:
This appication is an interactive pipeline in which users can leverage their domain knowledge to dynamically fix errors from the mutation callers.
The pipeline will take DNA sequencing data, perform mutation calling and processing to obtain a list of mutations in each tumor, and then run an algorithm for identifying drivers.
Users will be able to view the raw and processed data, edit it to address false positives/negatives, and view the resulting predicted drivers.
Furthermore, because the datasets are large, various strategies were employed for quickly running the different algorithms. 

How to execute the application:
Sample input files can be found in data folder.
To run the webserver, go to web folder, and execute python webserver.py.
Then, point your browser to 127.0.0.0:8080 to run the client.
