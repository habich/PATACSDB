#!/usr/bin/python

from flask import Flask

app = Flask(__name__)
app.config["SQLALCHEMY_DATABASE_URI"] = 'sqlite:///AA.db'
app.config["APPLICATION_ROOT"] = '/patacsdb'

import PATACSDB.views

if __name__ == "__main__":
    app.run()
