#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 25 2:20 PM 2023
Created in PyCharm
Created as QGP_Scripts/job_finshed_alert_whatsapp

@author: Dylan Neff, Dylan
"""

import os
import time

from WAMessenger import WAMessenger


def main():
    watch_jobs()
    # test_messenger_hello_world()
    print('donzo')


def watch_jobs():
    check_period = 60  # s
    job_checker = RCFJobChecker()
    jobs_alive = job_checker.get_total_jobs()
    while jobs_alive > 0:
        print(f'{jobs_alive} jobs still alive. Waiting {check_period} seconds before checking again')
        time.sleep(check_period)
        jobs_alive = job_checker.get_total_jobs()

    messenger = WAMessenger()
    # res = messenger.send_message('No more jobs running')
    res = messenger.send_message(template_name='hello_world')
    if res.ok:
        print('Jobs finished and WhatsApp message sent')
    else:
        print(f'Message sending failed: {res.json()}')


def test_job_checker():
    checker = RCFJobChecker()
    print(checker.get_total_jobs())


def test_messenger():
    messenger = WAMessenger()
    res = messenger.send_message('Test message')
    print(res.json())


def test_messenger_hello_world():
    messenger = WAMessenger()
    res = messenger.send_message(template_name='hello_world')
    print(res.json())


class RCFJobChecker:
    def __init__(self):
        self.user = 'dneff'

    def get_total_jobs(self):
        response = os.popen(f'condor_q {self.user} | tail -3')
        total_jobs = int(response.read().split('\n')[0].split(' jobs;')[0].replace('Total for query: ', ''))

        return total_jobs




# def get_phone_number(phone_number_path):
#     with open(phone_number_path, 'r') as file:
#         key, phone_number = file.readlines()
#     cipher_suite = Fernet(bytes(key, 'utf-8'))
#
#     return cipher_suite.decrypt(bytes(phone_number, 'utf-8')).decode('utf-8')


if __name__ == '__main__':
    main()
