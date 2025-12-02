import os
import shutil
import gdown

# --- CONFIGURATION ---------------------------------------------------------
# The Google Drive Folder URL
DRIVE_FOLDER_URL = 'https://drive.google.com/drive/folders/1UIN-r0rtcU1xaMbgnn450Pr3LcsqaKFs?usp=sharing'

# REPLACE THESE PLACEHOLDERS with the actual names of the files inside the Drive folder
# Ensure you keep the file extensions (e.g., 'data_2023.csv')
files_to_raw    = ["file1.csv", "file2.csv"]        # These go to /_raw
files_to_data   = ["file3.csv"]                     # These go to /data
files_to_output = ["file4.csv", "file5.csv"]        # These go to /data/output

# Temporary folder to hold the download before sorting
TEMP_DIR = "temp_drive_download"
# ---------------------------------------------------------------------------

def create_dirs():
    """Create the required directory structure."""
    dirs = ["_raw", "data", "data/output"]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
        print(f"✅ Directory ready: {d}")

def download_drive_folder():
    """Download the Google Drive folder contents."""
    print(f"\n⬇️  Downloading folder from Google Drive...")
    
    # Check if temp dir exists and clean it if needed to avoid conflicts
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
    
    # Download the folder using gdown
    gdown.download_folder(url=DRIVE_FOLDER_URL, output=TEMP_DIR, quiet=False, use_cookies=False)
    print("✅ Download complete.")

def move_files():
    """Move files from temp folder to their specific destinations."""
    print("\nCc️  Moving files to project structure...")

    # Dictionary mapping destination folders to the list of files intended for them
    moves = {
        "_raw": files_to_raw,
        "data": files_to_data,
        "data/output": files_to_output
    }

    for dest_folder, file_list in moves.items():
        for filename in file_list:
            src_path = os.path.join(TEMP_DIR, filename)
            dst_path = os.path.join(dest_folder, filename)

            if os.path.exists(src_path):
                shutil.move(src_path, dst_path)
                print(f"moved: {filename} -> {dest_folder}/")
            else:
                print(f"⚠️  WARNING: Could not find '{filename}' in the downloaded folder. Check the name in CONFIGURATION.")

def cleanup():
    """Remove the temporary download folder."""
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
        print("\n🧹 Cleanup complete (temp folder removed).")

def main():
    create_dirs()
    download_drive_folder()
    move_files()
    cleanup()
    print("\n🚀 CLARITY data import finished successfully.")

if __name__ == "__main__":
    main()