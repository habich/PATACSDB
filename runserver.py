#!/usr/bin/python
from PATACSDB import app
from flask.ext.sqlalchemy import get_debug_queries

@app.after_request
def after_request(response):
    for query in get_debug_queries():
        app.logger.warning("QUERY: %s\nParameters: %s\nDuration: %fs\nContext: %s\n" % (query.statement, query.parameters, query.duration, query.context))
    return response

if __name__ == "__main__":
    app.debug=True
    app.run(host='0.0.0.0')
