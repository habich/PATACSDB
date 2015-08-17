Install: SQLAlchemy, Flask, Flask-SQLAlchemy 

Run download.sh to download data from ensemble database. It will destroy existing data (previous versions) and download current version (even if there is no change). 
Run parse_data.py to analyse existing data about organisms, create files with informations about founded polyA (not required for database) and to initialize database.
! Warning: parse_data.py destroys database content.
Run runserver.py to run the server.


