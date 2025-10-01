import synapseclient
import pandas as pd
import argparse
import os
import sys

def download_files_from_synapse(auth_token_filepath, id_filepath, download_location):
    """
    Logs into Synapse and downloads a list of files.

    Args:
        auth_token_filepath (str): Path to a file containing the Synapse auth token.
        id_filepath (str): The path to a file containing Synapse IDs (one per line).
        download_location (str): The local directory to save downloaded files.
    """
    # --- 1. Read Auth Token from File ---
    try:
        with open(auth_token_filepath, 'r') as f:
            auth_token = f.readline().strip()
        if not auth_token:
            print(f"❌ Error: The token file '{auth_token_filepath}' appears to be empty.")
            sys.exit(1)
    except FileNotFoundError:
        print(f"❌ Error: The authentication token file '{auth_token_filepath}' was not found.")
        sys.exit(1)

    # --- 2. Log into Synapse ---
    syn = synapseclient.Synapse()
    try:
        print("Logging into Synapse...")
        syn.login(authToken=auth_token)
        print("✅ Login successful!")
    except Exception as e:
        print(f"❌ Error: Failed to log into Synapse. Please check your auth token.\nDetails: {e}")
        sys.exit(1)

    # --- 3. Create download directory if it doesn't exist ---
    if not os.path.exists(download_location):
        print(f"Creating download directory at: {download_location}")
        os.makedirs(download_location)

    # --- 4. Read Synapse IDs from the specified file ---
    try:
        print(f"Reading Synapse IDs from: {id_filepath}")
        syn_ids_df = pd.read_csv(filepath_or_buffer=id_filepath, header=None)
    except FileNotFoundError:
        print(f"❌ Error: The file '{id_filepath}' was not found.")
        sys.exit(1)

    # --- 5. Loop through IDs and download files ---
    print("\nStarting file download process...")
    # Iterate over each ID in the first column of the dataframe
    for syn_id in syn_ids_df[0]:
        try:
            print(f"  -> Downloading: {syn_id}")
            syn.get(syn_id, downloadLocation=download_location)
        except Exception as e:
            # Print a warning but continue with the next file
            print(f"  ⚠️  Warning: Could not download {syn_id}. Reason: {e}")
            
    print("\n✨ All downloads attempted. Script finished!")

if __name__ == '__main__':
    # --- Set up the command-line argument parser ---
    parser = argparse.ArgumentParser(
        description="Download a list of files from Synapse programmatically."
    )
    
    parser.add_argument(
        '-a', '--PATFile',
        required=True,
        help='Path to a text file containing your Synapse personal access token.'
    )
    
    parser.add_argument(
        '-f', '--SynapseIDFile',
        required=True,
        help='Path to the .txt file containing Synapse IDs, one ID per line.'
    )
    
    parser.add_argument(
        '-d', '--DownloadLocation',
        required=True,
        help='The local folder path to download the files into.'
    )

    # Parse the arguments provided by the user
    args = parser.parse_args()
    
    # Call the main function with the user-provided arguments
    download_files_from_synapse(args.authTokenFile, args.filepath, args.downloadLocation)
