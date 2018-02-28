#!/usr/bin/env python3
"""
:mod:`email_sender` -- 
========================================

.. module email_sender
   :platform: Unix, Windows, Mac, Linux
   :synopsis: 
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from lazy_property import LazyWritableProperty


class EmailSender:
    def __init__(self, from_address: str, to_address: str, subject: str = '(No Subject)'):
        self.from_address = from_address
        self.to_address = to_address
        self.subject = subject

    @LazyWritableProperty
    def server(self):
        pass

    @LazyWritableProperty
    def from_address_password(self) -> str:
        pass

    @LazyWritableProperty
    def message(self):
        pass

    def send(self):
        message = MIMEMultipart(self.message)
        message['Subject'] = self.subject
        message['From'] = self.from_address
        message['To'] = self.to_address
        message.attach(MIMEText(_text=self.message, _subtype='plain'))
        server = smtplib.SMTP_SSL('smtp.gmail.com')
        server.login(self.from_address, self.from_address_password)
        server.sendmail(from_addr=self.from_address, to_addrs=self.to_address, msg=message.as_string())
        server.quit()
