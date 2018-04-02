#!/bin/bash
NAME=cdm
HOMEDIR=/home/sasha_alesin/Documents/bmstu/bushuev/conjdirmethod
DJANGODIR=${HOMEDIR}/${NAME}
NUM_WORKERS=3
DJANGO_WSGI_MODULE=cdm/${NAME}.wsgi
GUNICORN=${HOMEDIR}/cdm/venv/bin/gunicorn

cd $HOMEDIR
source cdm/venv/bin/activate

exec ${GUNICORN} ${DJANGO_WSGI_MODULE}:application \
  --workers $NUM_WORKERS \
  --bind=unix:$HOMEDIR/run/gunicorn.sock \
  --log-file ${HOMEDIR}/logs/gunicorn.log