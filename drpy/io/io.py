from __future__ import absolute_import
import sys
import subprocess
import os
import numpy as np
import datetime
import re
from copy import deepcopy
from urllib.parse import urljoin, urlparse

# Canonical PPS directory for GPM DPR Level-2 NRT (listing + granule URLs)
NRT_DPR_L2_BASE_URL = 'https://jsimpsonhttps.pps.eosdis.nasa.gov/radar/DprL2/'

# NRT directory layout and access patterns are described in (login required):
# https://jsimpsonhttps.pps.eosdis.nasa.gov/documentation/nrtInstructions.pdf
# Public overview: https://gpm.nasa.gov/data/sources/pps-nrt
#
# Research: ...V9-20211125.YYYYMMDD-S...-E....032332.V07A.HDF5
# Research: ...V10-20260310.YYYYMMDD-S...-E....068975.V08A.nc
# NRT:      ...V1020260310.YYYYMMDD-S...-E....V08A.RT-NC  (optional orbit segment before .VxxA)
_DPR_FILENAME_RE = re.compile(
    r'2A\.GPM\.DPR\.V\d+[^.]*\.(?P<date>\d{8})-S(?P<start>\d{6})-E(?P<end>\d{6})'
    r'(?:\.\d+)?\.[^.]+\.(?P<ext>HDF5|RT-H5|RT-NC|nc)$',
    re.IGNORECASE
)


def _clean_listing_token(token):
    token = token.strip().strip('\'"<>')
    if token.startswith('href='):
        token = token.split('=', 1)[1].strip('\'"')
    token = token.split('?', 1)[0].strip('\'"<>')
    return token


def _extract_file_list_from_stdout(stdout):
    href_matches = re.findall(r'href=["\']([^"\']+)["\']', stdout, flags=re.IGNORECASE)
    candidates = href_matches + stdout.split()
    file_list = []
    for candidate in candidates:
        clean = _clean_listing_token(candidate)
        basename = os.path.basename(clean)
        if _DPR_FILENAME_RE.search(basename):
            file_list.append(clean)
    return list(dict.fromkeys(file_list))


def _curl_listing(url, username):
    # -f: fail on HTTP errors (avoid treating an HTML error page as a file listing)
    # -L: follow redirects; -sS: silent but show errors
    args = [
        'curl', '-sS', '-f', '-L', '-u', username + ':' + username, url,
    ]
    cmd = ' '.join(args)
    process = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    return cmd, stdout.decode(), stderr.decode(), process.returncode


def _parse_dpr_times(file_entry):
    basename = os.path.basename(file_entry)
    match = _DPR_FILENAME_RE.search(basename)
    if match is None:
        return None
    date_str = match.group('date')
    s_time = datetime.datetime.strptime(date_str + match.group('start'), '%Y%m%d%H%M%S')
    e_time = datetime.datetime.strptime(date_str + match.group('end'), '%Y%m%d%H%M%S')
    if s_time > e_time:
        e_time = e_time + datetime.timedelta(days=1)
    return s_time, e_time


_HDF5_FILE_MAGIC = (b'\x89HDF\r\n\x1a\n', b'\x89HDF\n\r\x1a\n')


def _verify_downloaded_product(path, min_bytes=4096):
    """
    Reject tiny or HTML-looking payloads (common when auth fails but curl still wrote a body).
    Real DPR granules are much larger than a few KB.
    """
    if not os.path.isfile(path):
        raise RuntimeError('Download did not create file: {!r}'.format(path))
    size = os.path.getsize(path)
    with open(path, 'rb') as f:
        head = f.read(512)
    stripped = head.lstrip()
    name = os.path.basename(path).lower()
    if size < min_bytes:
        raise RuntimeError(
            'Downloaded file is only {} bytes; expected a large binary granule. '
            'Often this means HTTP auth or path failed (see PPS NRT registration and '
            'https://jsimpsonhttps.pps.eosdis.nasa.gov/documentation/nrtInstructions.pdf ).'
            .format(size)
        )
    if stripped.startswith(b'<!DOCTYPE') or stripped.startswith(b'<html') or stripped.startswith(b'<?xml'):
        raise RuntimeError(
            'Downloaded content looks like HTML/XML, not radar data. '
            'Confirm PPS login (email as password), NRT access on your account, and the file URL.'
        )
    if name.endswith(('.rt-nc', '.rt-h5', '.hdf5', '.hdf')):
        if not any(head.startswith(m) for m in _HDF5_FILE_MAGIC):
            raise RuntimeError(
                'Downloaded file does not have an HDF5 signature (expected for NRT .RT-NC / .RT-H5). '
                'The file may be corrupt, truncated, or not a real PPS product. '
                'See https://jsimpsonhttps.pps.eosdis.nasa.gov/documentation/nrtInstructions.pdf'
            )


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
    """Query and download GPM-DPR from PPS (Research or NRT).

    NRT server layout and conventions are documented in NASA's
    `nrtInstructions.pdf <https://jsimpsonhttps.pps.eosdis.nasa.gov/documentation/nrtInstructions.pdf>`_
    (requires browser login with your registered PPS email).
    """
    
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
            # Canonical NRT DPR L2 base (see NRT_DPR_L2_BASE_URL); /text mirror as fallback
            search_urls = [
                NRT_DPR_L2_BASE_URL,
                'https://jsimpsonhttps.pps.eosdis.nasa.gov/text/radar/DprL2/',
            ]
            if self.s_time is not None:
                year = padder(self.s_time.year)
                month = padder(self.s_time.month)
                day = padder(self.s_time.day)
                dir_str = year + '/' + month + '/' + day + '/radar/'
                search_urls.extend([
                    'https://jsimpsonhttps.pps.eosdis.nasa.gov/text/gpmdata/' + dir_str,
                    'https://jsimpsonhttps.pps.eosdis.nasa.gov/gpmdata/' + dir_str,
                ])
            file_list = []
            for url in search_urls:
                cmd, stdout, _stderr, _rc = _curl_listing(url, self.username)
                self.cmd = cmd
                file_list = _extract_file_list_from_stdout(stdout)
                if file_list:
                    break

        elif self.servername=='Research':
            year = padder(self.s_time.year)
            month = padder(self.s_time.month)
            day = padder(self.s_time.day)
            dir_str = year + '/' + month + '/' + day + '/radar/'
            search_urls = [
                'https://arthurhouhttps.pps.eosdis.nasa.gov/text/gpmdata/' + dir_str,
                'https://arthurhouhttps.pps.eosdis.nasa.gov/gpmdata/' + dir_str,
            ]
            file_list = []
            for url in search_urls:
                cmd, stdout, _stderr, _rc = _curl_listing(url, self.username)
                self.cmd = cmd
                file_list = _extract_file_list_from_stdout(stdout)
                if file_list:
                    break
                
        if not file_list:
            raise RuntimeError(
                'No DPR files found from listing request. '
                'Check PPS credentials, server path availability, or file naming updates.'
            )
        self.file_list = file_list 
        
        
    def get_file(self,username,filename,server='https://jsimpsonhttps.pps.eosdis.nasa.gov/text'):
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
        
        parsed_files = []
        dtimes_s = []
        dtimes_e = []
        for file_entry in self.file_list:
            parsed = _parse_dpr_times(file_entry)
            if parsed is None:
                continue
            parsed_files.append(file_entry)
            dtimes_s.append(parsed[0])
            dtimes_e.append(parsed[1])

        self.file_list = np.asarray(parsed_files,dtype=str)
        if self.file_list.size == 0:
            raise RuntimeError(
                'No parseable DPR filenames found in server listing. '
                'The server format may have changed.'
            )
        dtimes_s = np.asarray(dtimes_s,dtype='object')
        dtimes_e = np.asarray(dtimes_e,dtype='object')

        if (self.s_time is None) and (self.e_time is None):
            print('Warning, not time range selected. All filenames are returned')
            return self.file_list
        elif self.e_time is None:
            ind_l = np.where(dtimes_s <= self.s_time)
            ind_r = np.where(dtimes_e >= self.s_time)
            ind_b = np.intersect1d(ind_l,ind_r)
            self.filename = self.file_list[ind_b]
        elif self.servername=='NearRealTime':
            ind_l = np.where(dtimes_s >= self.s_time)
            ind_r = np.where(dtimes_e <= self.e_time)
            ind_b = np.intersect1d(ind_l,ind_r)
            self.filename = self.file_list[ind_b]
        else:
            ind_l = np.where(dtimes_e >= self.s_time)
            ind_r = np.where(dtimes_s <= self.e_time)
            ind_b = np.intersect1d(ind_l,ind_r)
            self.filename = self.file_list[ind_b]

        if len(self.filename) == 0:
            raise RuntimeError(
                'No DPR file found overlapping the requested time range.'
            )

    def _build_download_url(self,file):
        file = str(file)
        parsed = urlparse(file)
        if parsed.scheme in ('http','https'):
            return file

        host = urlparse(self.server)
        host_root = host.scheme + '://' + host.netloc + '/'
        file = file.lstrip('\'"')

        if file.startswith('/'):
            return urljoin(host_root,file.lstrip('/'))
        if file.startswith('gpmdata/') or file.startswith('radar/'):
            return urljoin(host_root, file)
        # NRT: listings often yield a bare filename; granules live under NRT_DPR_L2_BASE_URL
        # (not under /text/ — joining to self.server would produce .../text/2A.... which 404s)
        if self.servername == 'NearRealTime':
            return urljoin(NRT_DPR_L2_BASE_URL, os.path.basename(file))
        return urljoin(self.server.rstrip('/') + '/', file)

    def download(self,savedir='./'):
        if not hasattr(self, 'filename') or len(self.filename) == 0:
            raise RuntimeError('No files selected for download.')
        for i,file in enumerate(self.filename):
            url = self._build_download_url(file)
            if self.verbose:
                print('Downloading {} of {}: {}'.format(i+1,len(self.filename),url))

            out_path = os.path.normpath(
                os.path.join(savedir, os.path.basename(file))
            )
            parent = os.path.dirname(os.path.abspath(out_path))
            if parent and not os.path.isdir(parent):
                os.makedirs(parent, exist_ok=True)

            args = [
                'curl',
                '-sS',
                '-f',
                '-L',
                '-u',
                self.username + ':' + self.username,
                url,
                '-o',
                out_path,
            ]
            process = subprocess.Popen(
                args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=False,
                universal_newlines=True,
            )

            if self.verbose:
                it = -1
                while True:
                    output = process.stderr.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        if it == -1:
                            print(output.strip())
                            it = 0
                        else:
                            print(output.strip(), end=' \r')
                rc = process.poll()
            else:
                process.wait()
                rc = process.returncode

            if rc != 0:
                raise RuntimeError(
                    'curl download failed (exit {}). URL: {}\n'
                    'If the file on disk is small or HTML, check PPS credentials and NRT access; see '
                    'https://jsimpsonhttps.pps.eosdis.nasa.gov/documentation/nrtInstructions.pdf'
                    .format(rc, url)
                )
            _verify_downloaded_product(out_path)
        

        if self.verbose:
            print('Done')
