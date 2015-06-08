#!/usr/bin/python

from flask import Flask

app = Flask(__name__)
app.config["SQLALCHEMY_DATABASE_URI"] = 'sqlite:///AA.db'

import PATACSDB.views
