---
layout: default
---
# Week 2 lecture and tutorial

## <ins>**High-performance computing (HPC) clusters**<ins>

![server](/Images/Week02/server.png)

**Servers vs Clusters vs HPCs:** I will use these terms interchangably in this course. There are slight differences, but for our purposes we can treat them as the same.

![HPC](/Images/Week02/HPC_diagram.png)

## The Novus cluster:
The newest cluster at OSU in Corvallis is named Novus. It is administered by the Digital Research Infrastructure team.

[Help page for Novus](https://oregonstate.teamdynamix.com/TDClient/1935/Portal/KB/?CategoryID=6889)

System stats:
- 1,712 total cores
- 9.3TB of total memory
- 80TB of scratch storage (disk space).

<br />
<br />

HPCs allow us to run analyses that are highly parallel or require large amounts of memory and/or time to run.

![computing](/Images/Week02/comp.png)

Here is one example of how a very modest amount of parallelization improves the runtime of an analysis.

![ERCnet runtime](/Images/Week02/ERCnet.png)

<br />
<br />

## <ins>**Using SSH to log into the server**<ins>

**SSH (Secure Shell):**
- SSH is a cryptographic network protocol used for secure communication between a client (your computer) and a server (the HPC in Corvallis) over an unsecured network.
- In our common vernacular, we'll sometimes use SSH as a verb. For example, "SSH-ing into the HPC" means use the `ssh` command to access the Novus server. 

## Using the ssh command
To use SSH to access the server, run the following command:

```ssh <your ONID>@novus.dri.oregonstate.edu```

- `ssh` is the command
- `<your ONID>` is your ONID (this should be the same as whatever comes before @oregonstate.edu in your email address)
- `novus.dri.oregonstate.edu` is the IP address of the Novus server

The above command will ask for your password. Type your password and press [enter]. For security, it will not show the letters as you type in your password, which can be a little disorienting, but you'll get used to it!

Once you're SSH'd in, any commands you run in your terminal will be executed on the remote machine (the HPC). The HPC runs linux, so any of your linux commands (e.g. `ls`, `cd`, etc...) should work just fine. It can be helpful to have two terminal windows open, one for running commands on your local machine, and one for running commands on your remote machine. 

**Important note:** You can only log into Novus when you are on campus (and connected to campus wifi). If you need to SSH into Novus when you're off campus, you can connect via VPN (see the instructions for 'VPN-ing' onto the campus system).

## Using VPN for access off-campus

[Help page for using VPN to access campus resources (including the HPC) when off campus.](https://oregonstate.teamdynamix.com/TDClient/1935/Portal/KB/?CategoryID=6889)

- We won't cover installing the VPN service on your laptop, but I would recommend installing it so that you can work on homework from home.

## Transferring data to and from the HPC using the SCP (secure copy) command

### **To upload data to the HPC:**
To upload a file from your local machine to the hpc, run the following command (from the command line on your local machine):
```
scp <file-to-upload> <your ONID>@novus.dri.oregonstate.edu:<path-to-where-you-want-the-file-to-land>
```
- `scp` is the command
- `<file-to-upload>` is the file to upload. To run it this way, you need to cd to the directory where this file lives.
- `<your ONID>@novus.dri.oregonstate.edu` should be the same from you ssh command.
- `:` is needed to specify a location on the remote machine.
- `<path-to-where-you-want-the-file-to-land>` is the full path to an existing directory on the remote machine. This needs to be a directory where you have 'write-permissions'. If you don't include the path, the file will be uploaded to your home directory by default.

Here is an example of how I would upload a file from my local machine to a directory called "test_upload_dir":
```
scp test_file.txt forsythe@novus.dri.oregonstate.edu:~/test_upload_dir/
```
- Note that scp is similar to the `cp` command:
  - `test_file.txt` is the 'source'
  - `forsythe@novus.dri.oregonstate.edu:~/test_upload_dir/` is the 'destination'

### **To download data from the HPC:**
```
scp <your ONID>@novus.dri.oregonstate.edu:<path-to-the-file-you-want-to-download> .
```
- `scp` is the command
- `<your ONID>@novus.dri.oregonstate.edu` tells scp where to look for the remote file
- `:` is needed to specify a location on the remote machine
- `<path-to-the-file-you-want-to-download>` is the full path to an existing file on the remote machine. This path should **include the name of the file**
- `.` means 'put the file in my current working directory (on the local machine)'

Here is an example of how I would download a file (~/test_download_dir/test_download_file.txt) from the remote machine to my local machine. 
```
scp forsythe@novus.dri.oregonstate.edu:~/test_download_dir/test_download_file.txt .
```
- Note that now:
  - `forsythe@novus.dri.oregonstate.edu:~/test_download_dir/test_download_file.txt` is the 'source'
  - `.` is the 'destination'


## Compressing and decompressing files for faster upload/download speed

### Creating a tar.gz File

**Create the tar.gz file**: Use the `tar` command to create the tar.gz archive. The basic syntax is:

   ```bash
   tar -czvf <name-of-zipped-file-to-create>.tar.gz <name-of-file-or-dir-to-compress>
   ```

   - `-c`: Create a new archive.
   - `-z`: Compress the archive with gzip.
   - `-v`: Verbose mode to display the files being archived.
   - `-f`: Specifies the filename of the archive.

### Decompressing a tar.gz File

**Decompress the tar.gz file**: Use the `tar` command with the `-x` option to extract the contents of the tar.gz archive.

   ```bash
   tar -xzvf </path/to/your/zipped-file.tar.gz>
   ```

   - `-x`: Extract files from the archive.
   - `-z`: Decompress the archive with gzip.
   - `-v`: Verbose mode to display the files being extracted.
   - `-f`: Specifies the filename of the archive.


<br />
<br />
## <ins>**Navigating the HPC**<ins>

File storage on the HPC is similar to any other file storage system. You can navigate the file system with standard commands like, `cd`, `ls`, and `pwd`.

There are a few key differences that are specific to HPCs/clusters:
- **Home directory (~):** your home directory is located in: `/home/<ONID>`. This is where you can store your files long-term. It's also a good place to install software. <ins>IMPORTANT: your home directory has a very small storage allotment. It is not intended as a place to put large data files.<ins>
- **Scratch storage:** The HPC is equipped with a region of the disk that is intended for storage of large amounts of data. Each user has their own allotment of 'scratch storage', located in `/scratch/<ONID>`. We'll call this directory 'your scratch directory'. Storage space here is virtually unlimited. <ins> However, scratch storage is temporary; scratch directories are deleted every ~3 months.<ins>
- **Shared directory for this class:**: I have created a shared directory for anyone in our group. Note that you may only have 'read permissions' for files here. If you need to edit a file, you'll have to copy the file into your home or scratch directory. 

## File permissions on the HPC
In a Linux system, permissions are represented by a 10-character string. You can see the permissions of any file/directory with `ls -l <file-or-dir-name>`. 

The permission string for files in the shared directory are: `-rw-r--r--`.

Here's what each character signifies:

1. **File Type Indicator**: The first character indicates the type of file. In this case, since it's a hyphen (`-`), it denotes a regular file. Another possible values= for this position would be `d` for directory.

2. **Owner Permissions**: Characters 2-4 (`rw-`) represent permissions for the owner of the file. In this example, the owner has read (`r`) and write (`w`) permissions, but not execute (`x`) permissions.

3. **Group Permissions**: Characters 5-7 (`r--`) represent permissions for the group that the file belongs to. In this case, the group has only read (`r`) permission.

4. **Other Permissions**: Characters 8-10 (`r--`) represent permissions for all other users who are not the owner and are not part of the group. In this example, other users also have only read (`r`) permission.



<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 1: compress a sequence file and download to your local machine.</strong>
   <ol>
      <li>SSH into the HPC</li>
      <li>Locate the fasta file in the class shared directory.</li>
      <li>Create a compressed (tar.gz) version of the file</li>
      <li>Use scp to download the tar.gz version of the file to your local machine.</li>
   </ol>
</div>

<br />
<br />



## <ins>**Linux text editors**<ins>

It is important to be able to edit text files (e.g. python scripts) in linux. When accessing the HPC via `ssh`, we'll need to use a text edit that is not a GUI. There are many linux text editors (e.g. `vim`, `nano`, `emacs`). The choice of which one to use comes down to what's available in your system and your personal preference. Below is a brief into to `vim`, but you're welcome to use any editor you prefer.

## Vim text editor
Vim is a powerful text editor that operates from the command line. It has a steep learning curve but offers efficient text editing once you get the hang of it. This tutorial will cover some basic commands to get you started with Vim.

### Opening Vim

To open a file in Vim, simply type `vim` followed by the name of the file you want to edit:

```bash
vim filename.txt
```

### Modes

Vim operates in different modes, each serving a distinct purpose:

1. **Command Mode**: Default mode for navigation and executing commands.
2. **Insert Mode**: Mode for inserting and editing text.
3. **Visual Mode**: Mode for selecting blocks of text.

To switch between modes:

- **Command Mode**: Press `Esc`.
- **Insert Mode**: Press `i` while in Normal Mode.
- **Visual Mode**: Press `v` while in Normal Mode.

### Navigation

In Command Mode, you can navigate through the text using the following commands:

- **Arrow keys**: Move the cursor in the respective direction.
- **h**: Move left.
- **j**: Move down.
- **k**: Move up.
- **l**: Move right.
- **Ctrl+f**: Move forward one page.
- **Ctrl+b**: Move backward one page.
- **gg**: Move to the beginning of the file.
- **G**: Move to the end of the file.
- **:x**: Move to a specific line number (replace `x` with the line number).

### Editing

In Insert Mode, you can type and edit text directly into the file.

- **i**: Enter Insert Mode before the cursor.

### Saving and Quitting

To save changes and exit Vim:

1. Press `Esc` to ensure you are in Command Mode.
2. Type `:wq` and press `Enter`.

To quit Vim without saving changes:

1. Press `Esc` to ensure you are in Command Mode.
2. Type `:q!` and press `Enter`.

![computing](/Images/Week02/vim.png)

<br />
<br />

<div style="border: 1px solid black; padding: 10px; margin: 10px 0;">
   <strong>Task 2: use vim to create a python script</strong>
   <ol>
      <li>Use vim to create a new python script</li>
      <li>Edit the script to get it to print "Hello world!" and then "Hello, again!" 100 times.</li>
      <li>Make the script 'executable' and then run it using './your_script.py'. </li>
      <li>Store the output of your script to a file using 'redirect' (>). </li>
   </ol>
</div>

<br />
<br />




<br />
<br />
## <ins>**VScode access to the server**<ins>
Virtual Studio Code (VScode) is an application you ca install on your local machine to provide a GUI (ish) environment for navigating the command line, editing files, etc... There is an SSH plugin that allows you to use VScode to access the HPC!

## VScode tutorial:

1. **Install Visual Studio Code**: If you haven't already, download and install Visual Studio Code from the [official website](https://code.visualstudio.com/).

2. **Install Remote Development Extension Pack (Optional)**: To facilitate SSH connections and remote development, it's recommended to install the "Remote Development" extension pack. You can install it by going to the Extensions view in VSCode (`Ctrl+Shift+X`), searching for "Remote Development", and clicking on "Install".

3. **Open a Project or Folder**: Open Visual Studio Code and either create a new project or open an existing folder where your code files are located.

4. **Open Remote SSH**: Once the "Remote Development" extension pack is installed, you'll have access to the remote development features. Click on the green button in the bottom-left corner of the status bar that says "Open a Remote Window". Then select "Remote-SSH: Connect to Host..." from the dropdown menu.

5. **Add SSH Host**: If you haven't configured any SSH hosts before, you'll need to add one. Click on "Add New SSH Host..." and enter the SSH connection details for your remote server. This includes the hostname (or IP address), username, and optional SSH key file.

6. **Connect to Remote Server**: After adding the SSH host, it will appear in the list of available hosts. Click on the host you want to connect to. VSCode will establish an SSH connection to the remote server.

7. **Authenticate (if necessary)**: If this is your first time connecting to the server, you may be prompted to authenticate the SSH connection. Depending on your server configuration, you may need to enter a password or passphrase, or you might be asked to confirm the server's fingerprint.

8. **Access Remote Terminal**: Once connected, VSCode will open a new window with access to the remote server. You'll have access to a terminal within VSCode, allowing you to execute commands directly on the remote server.

9. **Use Remote Server**: You can now work with files on the remote server directly within VSCode. Any changes you make to files will be reflected on the remote server.


<br />
<br />
## <ins>**Installing software with conda**<ins>

<br />
<br />
## <ins>**Running large jobs with SLURM job scheduler**<ins>

<br />
<br />
## <ins>**Using git and Github for collaborative code development**<ins>





