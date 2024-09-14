### New Zealand COVID-19 vaccine data (based on data released by Barry Young) 
Analysis of New Zeland COVID-19 vaccine data leaked by Barry Young to a US journalist who then made it public.

The [data](data) directory has the data released by Barry Young (obfuscated to remove any PII and dates have been shifted slightly so no PII is released).

The [code](code) directory has the file to run to do the batch analysis. Just do `python vax.py` after you `cd` to the `code` directory.

The [analysis](analysis) directory has an excel file where we've used the output data in pivot tables to show the 1 year mortality rates for each batch for each age range.

### Summary of vax.py

ChatGPT offers this helpful summary of vax.py:

1. **Purpose**: This code analyzes a dataset of vaccine records. It checks if people died within a year after receiving a vaccine and summarizes the results.
    
2. **Steps in the Code**:
    
    - **Read the Data**: The code reads a file containing vaccine records into a table (DataFrame).
    - **Format Dates**: It changes the date columns to a format that the code can work with easily.
    - **Add a New Column**: It creates a new column to show the month and year when the vaccine was given.
    - **Calculate Time Difference**: It figures out how many days passed between the vaccine date and the date of death.
    - **Determine Deaths within a Year**: It adds a column that indicates if someone died within a year of receiving the vaccine.
    - **Group and Summarize Data**: It groups the data by age, vaccine batch, dose number, and month-year of service. It counts the number of records and sums up the number of deaths within a year for each group.
    - **Save the Results**: Finally, it saves this summarized information into a new CSV file.
3. **Outcome**: The result is a new file that shows how many people received each batch of the vaccine, how many doses they got, and how many of those people died within a year of getting the vaccine, organized by age and month-year.
    
In essence, the code is taking a big list of NZ vaccine records, checking who died within a year after vaccination, summarizing the results, and then saving this summary into a new file.

### Note
The data here is NOT the original NZ data, but it is a derivative work. 

It contains no PII that matches any living or dead person anywhere in the world.

It does not violate any US law including copyright law.