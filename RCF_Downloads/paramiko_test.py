#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 05 7:06 PM 2021
Created in PyCharm
Created as QGP_Scripts/paramiko_test

@author: Dylan Neff, Dyn04
"""


import os
import paramiko
import traceback
import socket
import sys
import getpass

hostname = 'rftpexp.rhic.bnl.gov'
username = 'dneff'
rsa_path = 'C:/Users/Dyn04/.ssh/id_rsa'
port = 22

def byte_count(xfer, to_be_xfer):
    print(f'transferred: {((xfer / to_be_xfer) * 100):.2f} %')

try:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect((hostname, port))
except Exception as e:
    print("*** Connect failed: " + str(e))
    traceback.print_exc()
    sys.exit(1)

try:
    t = paramiko.Transport(sock)
    try:
        t.start_client()
    except paramiko.SSHException:
        print("*** SSH negotiation failed.")
        sys.exit(1)

    try:
        keys = paramiko.util.load_host_keys(
            os.path.expanduser("~/.ssh/known_hosts")
        )
    except IOError:
        try:
            keys = paramiko.util.load_host_keys(
                os.path.expanduser("~/ssh/known_hosts")
            )
        except IOError:
            print("*** Unable to open host keys file")
            keys = {}

    # check server's host key -- this is important.
    key = t.get_remote_server_key()
    if hostname not in keys:
        print("*** WARNING: Unknown host key!")
    elif key.get_name() not in keys[hostname]:
        print("*** WARNING: Unknown host key!")
    elif keys[hostname][key.get_name()] != key:
        print("*** WARNING: Host key has changed!!!")
        sys.exit(1)
    else:
        print("*** Host key OK.")

    try:
        key = paramiko.RSAKey.from_private_key_file(rsa_path)
    except paramiko.PasswordRequiredException:
        password = getpass.getpass("RSA key password: ")
        key = paramiko.RSAKey.from_private_key_file(rsa_path, password)
    t.auth_publickey(username, key)

    if not t.is_authenticated():
        print("*** Authentication failed. :(")
        t.close()
        sys.exit(1)

    sftp = paramiko.SFTPClient.from_transport(t)

    print(sftp.listdir("."))

    sftp.get('/gpfs01/star/pwg/dneff/data/BES1/trees/output/11GeV/auau_B440BDC7A77ADD9F6D9898840304B9CE_199.root',
             'C:/Users/Dyn04/Desktop/auau11_test.root', callback=byte_count)
    # with sftp.open('/gpfs01/star/pwg/dneff/data/BES1/trees/output/11GeV/auau_B440BDC7A77ADD9F6D9898840304B9CE_199.root',
    #                'r') as f:
    #     data = f.read()

    t.close()

except Exception as e:
    print("*** Caught exception: " + str(e.__class__) + ": " + str(e))
    traceback.print_exc()
    try:
        t.close()
    except:
        pass
    sys.exit(1)