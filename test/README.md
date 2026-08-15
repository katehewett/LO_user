This is a test scenario to run a simple s5cmd while using (and then not using) cron

copy_file.py  ... is a program that will:   
* move test_file.txt from laptop (or klone) to kopah bucket liveocean-kmhewett   
* it also creates a logfile called s5_transfer.log  

When using cron, a second logfile, copy_file_run.log will be created.  
* cron jobs listed in laptop_test.txt and klone_test.txt  


NOTE: the fix is already in copy_file.py.  
But originally, I could not use s5cmd while using cron (see fix at lines 50-55, and read below). But, I could run the script when executing the python script in terminal without using crontabs. 

No error code appeared when trying to run using line command e.g. ```python copy_file.py```, or ```run copy_file.py```. However, when trying to use crontabs the file would not run, and I would get exit code 127 / File not found error is raised. 
And insdie s5_transfer.log I would see:

```
CRITICAL ERROR: 's5cmd' executable was not found in your system PATH.
s5cmd Exit Code: 127
```

I think -- when using cron -- because Python was trying to execute the bare command "s5cmd", but the directory where s5cmd was installed isn't included in the system PATH variable visible to Python at runtime.

This might have happened for two reasons:

1) Cron's Restricted Environment: When Python runs via cron, cron executes with a minimal default PATH (typically just /usr/bin:/bin). It doesn't automatically load your interactive terminal's PATH where tools like s5cmd usually live (e.g., /usr/local/bin, /opt/homebrew/bin, or ~/bin).

2) Subprocess Lookup Failure: subprocess.run(["s5cmd", ...]) relies on os.environ["PATH"] to locate binaries. Because s5cmd isn't in those standard paths, Python fails to launch the process and triggers the FileNotFoundError block in your try/except statement.

The Fix (see lines 50 - 55)
The original version of copy_file.py had line 51:
```
cmd = ["s5cmd", "cp", SOURCE_FILE, DEST_BUCKET]
```
This is when the error was encountered.  
To fix, I commented out line 51, and wrote lines 52-55:  
```
    if str(HOME_DIR) == '/Users/katehewett':
        cmd = ["/Users/katehewett/Applications/miniconda3/envs/loenv/bin/s5cmd", "cp", SOURCE_FILE, DEST_BUCKET]
    elif (str(HOME_DIR) == '/mmfs1/home/kmhewett'):
        cmd = ["/usr/local/bin/s5cmd", "cp", SOURCE_FILE, DEST_BUCKET]
```
This is a little slopy, but it works -- the paths are what came from typing 
``` which -a s5cmd ```
on my laptop (or klone)  

A better option might have bee: 
```
s5cmd_bin = shutil.which('s5cmd') or '/usr/local/bin/s5cmd'
cmd = [s5cmd_bin, "cp", SOURCE_FILE, DEST_BUCKET]
```
which is what I suggested to write in K_driver_roms00 and driver_roms00. 

This now works on my laptop and on klone while using cron.   

Maybe something happened during maintanence that may have 'upset' how the s5cmd is located on klone? Perhaps a custom system-wide /etc/environment override occured on klone, which left cron with its bare-bones default: PATH=/usr/bin:/bin.

A potential fix might include implementing:
s5cmd_bin = shutil.which('s5cmd') or '/usr/local/bin/s5cmd'
and then cmd = [s5cmd_bin,...  ]
as suggested in driver_roms00 K_driver_roms00 (or equivalent)  might fix our issue. 