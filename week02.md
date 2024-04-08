---
layout: default
---
# Week 1 lecture and tutorial

## <ins>**High-performance computing (HPC) clusters**<ins

**Servers vs Clusters vs HPCs:** I will use these terms interchangably in this course. There are slight differences, but for our purposes we can treat them as the same.

## The Novus cluster:
- 1,712 total cores
- 9.3TB of total memory
- 80TB of scratch storage (disk space).
- 


<br />
<br />
## <ins>**Using SSH to log into the server**<ins

**SSH (Secure Shell):** 
- SSHis a cryptographic network protocol used for secure communication between a client (your computer) and a server (the HPC in Corvallis) over an unsecured network.
- In our common vernacular, we'll sometimes use SSH as a verb. For example, "SSH-ing into the HPC" means use the `ssh` command to access the Novus server. 

## Using the ssh command
To use SSH to access the server, run the following command:
```ssh <your ONID>@novus.dri.oregonstate.edu```
-`ssh` is the command
-`<your ONID>` is your ONID (this should be the same as whatever comes before @oregonstate.edu in your email address)
- `novus.dri.oregonstate.edu` is the IP address of the Novus server

The above command will ask for your password. Type your password and press [enter]. For security, it will not show the letters as you type in your password, which can be a little disorienting, but you'll get used to it!

**Important note:** You can only log into Novus when you are on campus (and connected to campus wifi). For SSH-ing into Novus when you're off campus, see the instructions for 'VPN-ing' onto the campus system.

## Using VPN for access off-campus

[Help page for using VPN to access campus resources (including the HPC) when off campus.](https://oregonstate.teamdynamix.com/TDClient/1935/Portal/KB/?CategoryID=6889)

- We won't cover installing the VPN service on your laptop, but I would recommend installing it so that you can work on homework from home.

## Transferring data to and from the HPC using the SCP (secure copy) command

### **To upload data to the HPC: **
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

<br />
<br />
## <ins>**Navigating the server and accessing shared data**<ins>

<br />
<br />
## <ins>**Linux text editors**<ins>

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





