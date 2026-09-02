# README: github notes on LO forking and branching
This readme is a work in progress...  

It outlines the steps that Kate took to seperate + use Parker's LO as he retires.  

Goals: to turn off Parker's crontabs // keep LO organized and flexible // allow me to run thre forecast from LO (Kate's version) 

I = Kate; he = Parker.  

This README is part 1 of <span style="color:red">X</span> and includes details on how I forked and maintain my LO branch.

# Forking LO: 
Historically, I have forked Parker's LO repo and then submitted pull requests to make a change, and operated in my LO_user. We want to keep this setup for other users and me too... Adding to the setup is changing my forked LO repo, so that I can run the forecast. I will be making changes to LO. I also want to be able to submit pull requests to Parker's upstream LO version and have a seperate version of LO "Kate's Version"

Github rules:  
* I cannot create two direct forks of the same repository within a single personal GitHub account.  GitHub restricts each account to one personal fork per upstream repository.  

## Options of common workflows include 

### <span style="color:blue"> 1. Use Separate Branches *(My choice)* **  </span>
I will use separate branches inside my existing single fork to LO. And will:  

* Keep my fork's main branch synced with the original repository (Parker's LO).
* Create LO-branch-1 (and could even have LO-branch-2) from that main branch to keep my work isolated (so I can run the forecast).
* This allows me to easily open distinct pull requests to the upstream repo from each branch.

### 2. Fork into a GitHub Organization *(Maybe future and needs to be explored more)*
In the future, we can look into doing this. But not right now. There is an option to create a completly separate repository structure on GitHub, where I can create a free GitHub Organization. Steps may include:

* Creating a new GitHub Organization where Kate or Parker's personal account is the owner. 
* Navigate to the master repository and click Fork.
* Change the Owner dropdown from your personal username to your newly created organization.

### 3. Duplicate via the "Import Repository" Feature *(Not for Kate v.)*
If I want a second copy of the code to play around with, but I wasn't intending on to submitting pull requests back to the original author (here Parker), then I *could* duplicate LO. But this isn't a step I want to take here.  

*However, this is likely what other people do that are running some version of LO that they then edit and run, becasue they aren't showing up as 'forked' in Parker's LO repo.*

To do this, likely involves the following steps:

* While logged into [your] github.com, click the + icon in the top right of GitHub and choose Import repository.
* Paste the URL of the master repository.
* Give it a completely new name.
* This creates a normal, independent repository in your personal account. Note: Because GitHub treats this as a brand new repo rather than a true fork, you cannot use the built-in interface to open Pull Requests back to the original project. 

<span style="color:red">*This might be a good option for Parker to test Kate's LO version on klone? He could follow the steps below and execute Kate's LO version on klone adn see if it breaks... ***but I think it might make more sense to try removing his LO (master) on his klone's gscratch/macc/parker, and then clone Kate's forked LO repo to his /macc/parker folder on klone***... Or maybe name a branch that he can go back to then accept all my pull requests to merge Kates v with Parkers???* ***More on this later and after my LO version produces a working forecast.***</span>

Some sources:  
[Forks: https://docs.github.com](https://docs.github.com/en/pull-requests/reference/forks)  
[Forks: https://www.reddit.com](https://www.reddit.com/r/git/comments/y0hyl3/git_forkcontinually_changing_main_branch_and/)  
[Forks: https://www.reddit.com](https://www.reddit.com/r/github/comments/qvknv4/how_do_i_fork_a_repo_again/)   
[Forks: https://docs.github.com](https://docs.github.com/articles/fork-a-repo)    
[Double fork? Nope: https://stackoverflow.com](https://stackoverflow.com/questions/19279490/can-i-fork-another-persons-repo-twice-into-my-own-account)  
[ORGs: https://github.com](https://github.com/orgs/community/discussions/23095)   
[More ORGs info: https://github.com](https://docs.github.com/en/organizations/collaborating-with-groups-in-organizations/about-organizations)

## <span style="color:blue"> Option #1: Setup Separate Branches **</span>
On my laptop, I will track two copies of the code under different names:  
* upstream (the original LO repo that I only pull updates from) This is Parker's LO
* origin (my fork that I push edits to).  

Here is the exact step-by-step setup to achieve this:

### Step0. Fork the LO repo 
Create a fork of Parker's LO repo. 
To do this, I went to https://github.com/parkermac/LO
Clicked fork and followed directions

### Step1. Clone my forked LO repo 
Download my personal fork to my computer [and klone]
git clone https://github.com/katehewett/LO.git

This automatically set my fork as origin. 
When I run code or make changes, this is the copy that I edit.

To connect to the original LO repo, I 
copied the URL of the original parent repository from GitHub
https://github.com/parkermac/LO

### Step2. Connect to upstream repo

and then in terminal I did: 
```
git remote add upstream https://github.com/parkermac/LO
```

### Step3. Check branches and connections
In terminal, if I do $ (I see *):  
On laptop:
```
$ git branch -a
* main
  remotes/origin/HEAD -> origin/main
  remotes/origin/main
  remotes/upstream/HEAD -> upstream/main
  remotes/upstream/main
```

On klone:
```
 git branch -a
* main
  remotes/origin/HEAD -> origin/main
  remotes/origin/main
```
*Where Parker's LO = upstream*  
*my LO = origin*

And because I had connected the original LO repo, on my laptop when I do $ (I see): 
```
$ git remote -v
origin	https://github.com/katehewett/LO.git (fetch)
origin	https://github.com/katehewett/LO.git (push)
upstream	https://github.com/parkermac/LO.git (fetch)
upstream	https://github.com/parkermac/LO.git (push)
```

This output shows origin pointing to my fork (Kate's version), and upstream pointing to the original LO repo (Parker's version).

But on Klone, when I do $ (I see):
```
git remote -v
origin	https://github.com/katehewett/LO.git (fetch)
origin	https://github.com/katehewett/LO.git (push)
```

And that's because I had not connected to the upstream on Klone. <span style="color:red">(DO I WANT TO??? hold off for now)</span>

## <span style="color:blue"> How to's --  Seperate LO branches **</span>
Note, on my laptop:  

When I ```git remote -v``` in terminal, and I see upstream listed in my terminal output, that means my cloned repository is successfully connected to the original parent project (here Parker's LO).

Because Git handles the underlying configuration, GitHub Desktop will automatically detect this upstream connection the next time I open or refresh the app.  
If I open Github desktop / select my LO repository / and click on the Current Branch dropdown menu at the top center of the screen... then I will see that I have a Default Branch (currently set to) main, which is Kate's edited version of LO. And then under other branches (currently) see the upstream/main branch which is Parker's master LO.  

The cherry pick command will let me submit single pull requests to parker. 

can toggle back / forth and set up other named branches. but for now. staying put with these 2.