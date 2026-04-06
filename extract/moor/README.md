# README 
Contents:
* job list for multi mooring driver  
* dye_01 additions, see below

## dye_01 additions: 

These files were created and modified to help A Leeson with her dye text extractions for CREST. 
 
* dye01 can be included as a variable for your mooring extractions (see step1)
* I submitted a pull request in LO to add a list type -hisX which will let you pull the 2nd history files in your dye test (see step2)

#### Step1: These are found under Kate's LO_user:  
* extract_moor_dye_01.py  
* multi_mooring_driver_dye01.py
* jobs_list.py: see and add 'whidbey_test' if you want to replicte my comments section run 

You can add to your LO_user and it should run on your desktop and apogee. 

#### Step2: new list type, hisX

I added a new list type, under Lfun / get_fn_list. The default his_num is 2, but you can enter a different value if you don't want the defaulted 0002 history file. see below, new addition:

    elif list_type == 'hisX':
        # list of history files over a date range for history files associated with his_num entry. 
        # his_num default = 2, therefore will grab ocean_his_0002 over a date range of saved data
        his_string = ('0000' + str(his_num))[-4:]
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / ('ocean_his_' + his_string + '.nc') 
            #fn = dir0 / f_string / 'ocean_his_0002.nc'
            fn_list.append(fn)
















