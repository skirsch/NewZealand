'''
Analyze the NZ records file which is vaccine record level data.

'''
import pandas as pd
import csv  # for the quoting option on output
from datetime import timedelta

def process_data(file_path, output_file):
    df = pd.read_csv(file_path)

    # Convert date columns to datetime
    df['date_time_of_service'] = pd.to_datetime(df['date_time_of_service'], format='%m-%d-%Y')
    df['date_of_death'] = pd.to_datetime(df['date_of_death'], format='%m-%d-%Y')

    # Create a new column for month-year of service
    df['month_year_of_service'] = df['date_time_of_service'].dt.to_period('M')

    # Calculate whether the person died within 365 days
    time_diff = (df['date_of_death'] - df['date_time_of_service']).dt.days

    # Create a boolean column where True indicates death within 365 days, and False otherwise
    df['death_within_365_days'] = (time_diff <= 365) & (time_diff >= 0)

    # Group by and aggregate
    grouped_data = (df
                    .groupby(['age', 'batch_id', 'dose_number', 'month_year_of_service'])
                    .agg(
                        count=('mrn', 'count'),
                        deaths_within_365_days=('death_within_365_days', 'sum')  # Sum the boolean values (True = 1, False = 0)
                    )
    )

    # Reset the index to make the grouping columns regular columns
    grouped_data = grouped_data.reset_index()

    # Output the result to CSV, including the grouping columns
    grouped_data.to_csv(output_file, index=False)

# Call the function
file_path = '../data/nz-record-level-data-4M-records.csv'
output_file = '../data/analysis_by_batch.csv'
result = process_data(file_path, output_file)
