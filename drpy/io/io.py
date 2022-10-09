from __future__ import absolute_import
import sys
import subprocess
import os
import numpy as np
import datetime
from copy import deepcopy

def padder(x):
    
    if x < 10:
        x = '0' + str(x)
    else:
        x = str(x)
        
    return x 

def find_keys(keys,substring):
  """keys is a list of strings; substring is a list of substrings"""
  for sub in substring:
    res = [i for i in keys if sub in i]
    keys = deepcopy(res)

  return res

class netrunner():
    """ This class will house all the functions needed to query the GPM FTP"""
    
    def __init__(self,servername='NearRealTime',username=None,start_time=None,end_time=None,
                    autorun=True,savedir='./',verbose=True):
        self.servername = servername
        if servername=='NearRealTime':
            self.server ='https://jsimpsonhttps.pps.eosdis.nasa.gov/text'
        elif servername=='Research':
            self.server = 'https://arthurhouhttps.pps.eosdis.nasa.gov/text'
        self.s_time = start_time
        self.e_time = end_time
        self.verbose = verbose 

        #check username input 
        if username is None:
            print('Please enter your PPS registered email as the username')
        else:
            self.username=username
        
        #check dates, multi-day not supported 
        if self.e_time is not None:
            if self.s_time.day != self.e_time.day:
                print('More than one day given as input!')
                print('Multi-day downloads are not currently supported')

        if autorun:
            #this will grab all the files on your day of interest
            self.get_file_list()
            #this will grab the one file that has the time you gave it 
            self.locate_file()
            #this will download it locally
            self.download(savedir=savedir)
        
    def get_file_list(self):
        """ 
        Code to grab files that are on the servers to let the user decide which files they want
        """
        if self.servername=='NearRealTime':
            url = self.server + '/radar/DprL2/' 
            cmd = 'curl -s -u ' + self.username+':'+self.username+' ' + url
            self.cmd = cmd
            args = cmd.split()
            process = subprocess.Popen(args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
            stdout = process.communicate()[0].decode()
            file_list = stdout.split()
            file_list = find_keys(file_list,['2A.GPM.DPR.V920211125'])
            
        elif self.servername=='Research':
            server = 'https://arthurhouhttps.pps.eosdis.nasa.gov/text/'
            year = padder(self.s_time.year)
            month = padder(self.s_time.month)
            day = padder(self.s_time.day)
            dir_str = 'gpmdata/' + year + '/' + month + '/' + day + '/radar/' 
            url = server + dir_str
            cmd = 'curl -s -u ' + self.username+':'+self.username+' ' + url
            self.cmd = cmd
            args = cmd.split()
            process = subprocess.Popen(args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
            stdout = process.communicate()[0].decode()
            file_list = stdout.split()
            file_list = find_keys(file_list,['2A.GPM.DPR.V9-20211125'])
                
        self.file_list = file_list 
        
        
    def get_file(username,filename,server='https://jsimpsonhttps.pps.eosdis.nasa.gov/text'):
        """ Some bit of code modified from here: 
        https://gpm.nasa.gov/sites/default/files/document_files/PPS-jsimpsonhttps_retrieval.pdf
        """

        #Note from dev.: this function is not used... RJC 10/09/22
        url = server + filename

        if self.verbose:
            print('Downloading: {}'.format(url))

        cmd = 'curl -s -u ' + username+':'+username+' ' + url + ' -o ' + \
        os.path.basename(filename)
        args = cmd.split()
        process = subprocess.Popen(args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
        process.wait()

        if self.verbose:
            print('Done')

    def locate_file(self):
        """ 
        This is a method that will grab the desired GPM-DPR files from the NRT server
        username: your email you signed up for the PPS server (string)
        start_time: datetime of when you want to start looking 
        end_time: datetime of when you end looking
        """ 
        
        if self.servername=='NearRealTime':
            s_times = []
            e_times = []
            date = []
            for i in self.file_list:
                t = i.split('S')[1][:14].split('-E')
                s_times.append(t[0])
                e_times.append(t[1])
                #this is for the old version RJC 06-Dec-2021
                #date.append(i.split('/radar/DprL2/2A.GPM.DPR.V820180723.')[1].split('-S')[0])
                #new version string start RJC 06-Dec-2021
                date.append(i.split('/radar/DprL2/2A.GPM.DPR.V920211125.')[1].split('-S')[0])
                
                res = [i + j for i, j in zip(date, s_times)] 
                dtimes_s = np.zeros(len(res),dtype='object')
                res2 = [i + j for i, j in zip(date, e_times)] 
                dtimes_e = np.zeros(len(res),dtype='object')
            for i,ii in enumerate(res):
                dtimes_s[i] = datetime.datetime.strptime(ii,'%Y%m%d%H%M%S')
                dtimes_e[i] = datetime.datetime.strptime(res2[i],'%Y%m%d%H%M%S')

                if dtimes_s[i] > dtimes_e[i]:
                    dtimes_e[i] = dtimes_e[i] + datetime.timedelta(days=1)


            self.file_list = np.asarray(self.file_list,dtype=str)
            if (self.s_time is None) and (self.e_time is None):
                print('Warning, not time range selected. All filenames are returned')
                return self.file_list
            elif self.e_time is None:
                ind_l = np.where(dtimes_s <= self.s_time)
                ind_r = np.where(dtimes_e >= self.s_time)
                ind_b = np.intersect1d(ind_l,ind_r)
                self.filename = self.file_list[ind_b]
            else:
                ind_l = np.where(dtimes_s >= self.s_time)
                ind_r = np.where(dtimes_e <= self.e_time)
                ind_b = np.intersect1d(ind_l,ind_r)
                self.filename = self.file_list[ind_b]
            
        elif self.servername=='Research':
            
            s_times = []
            e_times = []
            date = []
            for i in self.file_list:
                t = i.split('S')[1][:14].split('-E')
                s_times.append(t[0])
                e_times.append(t[1])
                #V7 fix 
                date.append(i.split('2A.GPM.DPR.V9-20211125.')[1].split('-S')[0])

                res = [i + j for i, j in zip(date, s_times)] 
                dtimes_s = np.zeros(len(res),dtype='object')
                res2 = [i + j for i, j in zip(date, e_times)] 
                dtimes_e = np.zeros(len(res),dtype='object')
            for i,ii in enumerate(res):
                dtimes_s[i] = datetime.datetime.strptime(ii,'%Y%m%d%H%M%S')
                dtimes_e[i] = datetime.datetime.strptime(res2[i],'%Y%m%d%H%M%S')

                if dtimes_s[i] > dtimes_e[i]:
                    dtimes_e[i] = dtimes_e[i] + datetime.timedelta(days=1)


            self.file_list = np.asarray(self.file_list,dtype=str)
            if (self.s_time is None) and (self.e_time is None):
                print('Warning, not time range selected. All filenames are returned')
                return self.file_list
            elif self.e_time is None:
                ind_l = np.where(dtimes_s <= self.s_time)
                ind_r = np.where(dtimes_e >= self.s_time)
                ind_b = np.intersect1d(ind_l,ind_r)
                self.filename = self.file_list[ind_b]
            else:
                ind_l = np.where(dtimes_e >= self.s_time)
                ind_r = np.where(dtimes_s <= self.e_time)
                ind_b = np.intersect1d(ind_l,ind_r)
                self.filename = self.file_list[ind_b]

    def download(self,savedir='./'):
        for i,file in enumerate(self.filename):
            url = self.server + file
            if self.verbose:
                print('Downloading {} of {}: {}'.format(i,len(self.filename),url))

            cmd = 'curl -u ' + self.username+':'+self.username+' ' + url + ' -o ' + \
            os.path.basename(savedir+file)
            args = cmd.split()
            process = subprocess.Popen(args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=False,
            universal_newlines=True)

            # outs, errs = process.communicate()
            # print(outs,errs)
            # process.wait()
            if self.verbose:
                it = -1 
                while True:

                    output = process.stderr.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        if it == -1:
                            #need to print header 
                            print(output.strip())
                            it = 0
                        else:
                            print(output.strip(),end=' \r',)

                rc = process.poll()
            else:
                process.wait()
        

        if self.verbose:
            print('Done')
