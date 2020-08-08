from __future__ import absolute_import
import sys
import subprocess
import os
import numpy as np
import datetime
from ftplib import FTP

def padder(x):
    
    if x < 10:
        x = '0' + str(x)
    else:
        x = str(x)
        
    return x 

class netrunner():
    """ This class will house all the functions needed to query the GPM FTP"""
    
    def __init__(self,servername='NearRealTime',username=None,start_time=None,end_time=None,Xradar=False):
        self.servername = servername
        if servername=='NearRealTime':
            self.server ='https://jsimpsonhttps.pps.eosdis.nasa.gov/text'
        elif servername=='Research':
            self.server = 'pps.gsfc.nasa.gov'
        self.s_time = start_time
        self.e_time = end_time
        if username is None:
            print('Please enter your PPS registered email as the username')
        else:
            self.username=username
        self.Xradar = Xradar
        
        
    def get_file_list(self):
        """ 
        Code to grab files that are on the servers to let the user decide which files they want
        """
        if self.servername=='NearRealTime':
            url = self.server + '/radar/DprL2/' 
            cmd = 'curl -s -u ' + self.username+':'+self.username+' ' + url
            args = cmd.split()
            process = subprocess.Popen(args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
            stdout = process.communicate()[0].decode()
            file_list = stdout.split()
        elif self.servername=='Research':
            ftp = FTP(self.server)
            ftp.login(self.username,self.username)
            year = padder(self.s_time.year)
            month = padder(self.s_time.month)
            day = padder(self.s_time.day)
            if self.Xradar:
                dir_str = 'gpmdata/' + year + '/' + month + '/' + day + '/Xradar/' 
            else:
                dir_str = 'gpmdata/' + year + '/' + month + '/' + day + '/radar/' 
            ftp.cwd(dir_str)
            file_list = ftp.nlst()
            file_list.sort()
            self.ftp = ftp
            
            if self.Xradar:
                file_list = [i for i in file_list if '2A.GPM.DPRX.V8' in i]
            else:
                file_list = [i for i in file_list if '2A.GPM.DPR.V8' in i]
            
                
        self.file_list = file_list 
        
        
    def get_file(username,filename,server='https://jsimpsonhttps.pps.eosdis.nasa.gov/text'):
        """ Some bit of code modified from here: 
        https://gpm.nasa.gov/sites/default/files/document_files/PPS-jsimpsonhttps_retrieval.pdf
        """
        url = server + filename
        print('Downloading: {}'.format(url))
        cmd = 'curl -s -u ' + username+':'+username+' ' + url + ' -o ' + \
        os.path.basename(filename)
        args = cmd.split()
        process = subprocess.Popen(args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
        process.wait()
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
                date.append(i.split('/radar/DprL2/2A.GPM.DPR.V820180723.')[1].split('-S')[0])

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
                if self.Xradar:
                    date.append(i.split('2A.GPM.DPRX.V8-20200326.')[1].split('-S')[0])
                else:
                    date.append(i.split('2A.GPM.DPR.V8-20180723.')[1].split('-S')[0])

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

    def download(self,savedir='./'):
        
        if self.servername=='NearRealTime':
            url = self.server + self.filename[0]
            print('Downloading: {}'.format(url))
            cmd = 'curl -s -u ' + self.username+':'+self.username+' ' + url + ' -o ' + \
            os.path.basename(savedir+self.filename[0])
            args = cmd.split()
            process = subprocess.Popen(args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
            process.wait()
            print('Done')
        elif self.servername=='Research':
            print('Downloading: {}'.format(self.filename[0]))
            with open(savedir+self.filename[0], 'wb') as fp:
                self.ftp.retrbinary('RETR '+self.filename[0], fp.write)
            print('Done')