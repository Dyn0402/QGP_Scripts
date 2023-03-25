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

from cryptography.fernet import Fernet
import json
import requests
from time import sleep


def main():
    # watch_jobs()
    test_messenger_hello_world()
    print('donzo')


def watch_jobs():
    check_period = 60  # s
    job_checker = RCFJobChecker()
    while job_checker.get_total_jobs() > 0:
        time.sleep(check_period)

    messenger = WAMessenger()
    res = messenger.send_message('No more jobs running')
    print(res.json())


def test_job_checker():
    checker = RCFJobChecker()
    print(checker.get_total_jobs())


def test_messenger():
    messenger = WAMessenger()
    res = messenger.send_message('Test message')
    print(res.json())


def test_messenger_hello_world():
    messenger = WAMessenger()
    res = messenger.send_hello_world()
    print(res.json())


class WAMessenger:
    def __init__(self):
        self.data_paths = [
            'C:/Users/Dylan/Desktop/wamessenger_data.txt',
            '/star/u/dneff/wamessenger_data.txt'
        ]
        self.data_path = self.try_paths()
        self.key, self.phone_receive, self.phone_send, self.token = None, None, None, None
        self.cipher_suite = None
        self.read_data()

    def read_data(self):
        with open(self.data_path, 'r') as file:
            lines = [bytes(x, 'utf-8') for x in file.readlines()]
        self.key = lines.pop(0)
        self.cipher_suite = Fernet(self.key)
        self.phone_receive, self.phone_send, self.token = \
            [self.cipher_suite.decrypt(x).decode('utf-8') for x in lines[:3]]

    def send_message(self, message='Hello from RCF'):
        url = f'https://graph.facebook.com/v16.0/{self.phone_send}/messages'
        headers = {
            'Authorization': f'Bearer {self.token}',
            'Content-Type': 'application/json'
        }

        msg_body_params = [
            {
                "type": "text",
                "text": message
            }
        ]

        data = {
            'messaging_product': 'whatsapp',
            'to': self.phone_receive,
            'type': 'template',
            'template': {
                'name': 'rcf_message',
                'language': {
                    'code': 'en_US'
                },
                'components': [
                    {
                        'type': 'body',
                        'parameters': msg_body_params
                    }
                ]
            }
        }

        return requests.post(url, headers=headers, data=json.dumps(data))

    def send_hello_world(self):
        url = f'https://graph.facebook.com/v16.0/{self.phone_send}/messages'
        headers = {
            'Authorization': f'Bearer {self.token}',
            'Content-Type': 'application/json'
        }

        data = {
            'messaging_product': 'whatsapp',
            'to': self.phone_receive,
            'type': 'template',
            'template': {
                'name': 'hello_world',
                'language': {
                    'code': 'en_US'
                }
            }
        }

        return requests.post(url, headers=headers, data=json.dumps(data))

    def try_paths(self):
        """
        Try all data file paths until one is successfully opened
        :return:
        """
        good_path = None
        for path in self.data_paths:
            try:
                open(path, 'r')
            except FileNotFoundError:
                continue
            good_path = path
            break
        return good_path


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
