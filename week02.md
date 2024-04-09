---
layout: default
---
# Week 1 lecture and tutorial

## <ins>**High-performance computing (HPC) clusters**<ins>

**Servers vs Clusters vs HPCs:** I will use these terms interchangably in this course. There are slight differences, but for our purposes we can treat them as the same.

![HPC](/Images/Week02/HPC_diagram.png)

## The Novus cluster:
- 1,712 total cores
- 9.3TB of total memory
- 80TB of scratch storage (disk space).

![server](/Images/Week02/server.png)

<br />
<br />

![computing](/Images/Week02/comp.png)

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

<br />
<br />
## <ins>**Navigating the HPC**<ins>

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

1. **Normal Mode**: Default mode for navigation and executing commands.
2. **Insert Mode**: Mode for inserting and editing text.
3. **Visual Mode**: Mode for selecting blocks of text.

To switch between modes:

- **Normal Mode**: Press `Esc`.
- **Insert Mode**: Press `i` while in Normal Mode.
- **Visual Mode**: Press `v` while in Normal Mode.

### Navigation

In Normal Mode, you can navigate through the text using the following commands:

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

1. Press `Esc` to ensure you are in Normal Mode.
2. Type `:wq` and press `Enter`.

To quit Vim without saving changes:

1. Press `Esc` to ensure you are in Normal Mode.
2. Type `:q!` and press `Enter`.

![computing](/Images/Week02/vim.png)

<br />
<br />
## <ins>**VScode access to the server**<ins>

<br />
<br />
## <ins>**Installing software with conda**<ins>

<br />
<br />
## <ins>**Running large jobs with SLURM job scheduler**<ins>

<br />
<br />
## <ins>**Using git and Github for collaborative code development**<ins>





