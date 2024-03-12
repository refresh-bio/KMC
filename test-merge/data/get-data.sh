#!/bin/bash

if [ -e data.tar.gz ]; then
	echo "already downloaded"
else
    wget -O data.tar.gz --no-check-certificate  "https://mkokot.ddns.net:5001/d/s/xCVEwWDC2L7vmVonajfvX5Mg7TuroVn8/webapi/entry.cgi/SYNO.SynologyDrive.Files/data.tar.gz?api=SYNO.SynologyDrive.Files&method=download&version=2&files=%5B%22id%3A801418206364350490%22%5D&force_download=true&sharing_token=%22gPGbHjZmq6HO4d6NOXGm3dDLcclW4LtxKMc0DqSoHhFFh8WOlTi43MjigLnih7h1oU8chQOKsF8WxeR2w1neXstfSkcBXBmkJwE53.A7dN0v218kJvVYHrytKgrqw1xCiuYowzv_AAHwMKs7dBNwG3A2E5fV5YWwiaRTuU_2CHeqiWxMe5TwU4GcTt6VNJJaTivpHCqDjIlQfJA_Joq0fyRxnEzpSCYxEuGpb64Dlg59cSW0f9CgEGkj%22"
	tar -xvf data.tar.gz
fi
