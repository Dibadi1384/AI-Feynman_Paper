AI Feynman Data Generation & Execution
Overview
This project generates example datasets, runs AI Feynman on them, and collects results with runtime statistics.

Usage Instructions
1. Generate Data
Create a folder named example_data:

bash
Copy
Edit
mkdir example_data
Run the data generation script:

bash
Copy
Edit
python Dataget.py
Move the generated data files into the example_data folder:

bash
Copy
Edit
mv <generated_files> example_data/
2. Run AI Feynman
Execute:

bash
Copy
Edit
./aifeynman_job.sh
This will:

Run AI Feynman for all data files listed in feynex.py

Print the runtime for each dataset at the end

3. Modify Data Files to Run
To change which datasets AI Feynman processes, edit:

bash
Copy
Edit
feynex.py
4. View Results
Results are automatically saved in the results/ folder.


