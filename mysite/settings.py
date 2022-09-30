"""
Django settings for mysite project.

Generated by 'django-admin startproject' using Django 4.0.1.

For more information on this file, see
https://docs.djangoproject.com/en/4.0/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/4.0/ref/settings/
"""
import os
import yaml
import logging
from pathlib import Path
import pymysql 
pymysql.install_as_MySQLdb()

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/4.0/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY_FILE = 'django_secret_key' if os.path.exists('django_secret_key') else os.environ.get('DJANGO_SECRET_KEY_FILE')
with open(SECRET_KEY_FILE, 'r') as dsk_fh:
    SECRET_KEY = dsk_fh.read()

# SECURITY WARNING: don't run with debug turned on in production!
if os.environ.get('DEBUG_MODE') is None:
    DEBUG = True
else:
    DEBUG = (os.environ.get('DEBUG_MODE') == "True")    # This environment variable is defined in the Spin container


# ALLOWED_HOSTS = []
ALLOWED_HOSTS = [
    'biomoleculararchitect.lbl.gov',
    'molinv.molinv.development.svc.spin.nersc.org',
    '127.0.0.1', 
    'localhost', 
]



# Application definition

INSTALLED_APPS = [
    'retroapp.apps.RetroappConfig',
    'renderer.apps.RendererConfig',
    'django.contrib.admin',
    'django.contrib.sites',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'widget_tweaks',
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'mysite.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'mysite.wsgi.application'


# Database
# https://docs.djangoproject.com/en/4.0/ref/settings/#databases
db_config_file = "mysql_config.yaml"
if os.path.exists(db_config_file):
    logging.info(f"Logging to MySQL using {db_config_file} file.")
    with open(db_config_file, "r") as dbconf_fh:
        db_config = yaml.safe_load(dbconf_fh)
else:
    logging.info(f"Logging to MySQL using 'MYSQL_PASSWORD_FILE' environment variable.")
    with open(os.environ.get('MYSQL_PASSWORD_FILE'), 'r') as mysql_pass_fh:
        mysql_password = mysql_pass_fh.read()
    
    db_config = dict(
        DBNAME=os.environ.get('MYSQL_DATABASE'),
        USER=os.environ.get('MYSQL_USER'),
        PASSWORD=mysql_password,
        HOST=os.environ.get('MYSQL_HOST'),
        PORT=os.environ.get('MYSQL_PORT'),
    )
DATABASES = {
    # 'default': {
    #     'ENGINE': 'django.db.backends.sqlite3',
    #     'NAME': BASE_DIR / 'db.sqlite3',
    # }
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': db_config['DBNAME'],
        'USER': db_config['USER'],
        'PASSWORD': db_config['PASSWORD'],
        'HOST': db_config['HOST'],
        'PORT': db_config['PORT'],
        'OPTIONS': {
            'init_command': "SET sql_mode='STRICT_TRANS_TABLES'"
         }   
    }
}


# Password validation
# https://docs.djangoproject.com/en/4.0/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/4.0/topics/i18n/

LANGUAGE_CODE = 'en-us'

# TIME_ZONE = 'UTC'
TIME_ZONE = 'America/Los_Angeles'

USE_I18N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/4.0/howto/static-files/

STATIC_URL = 'static/'

# Default primary key field type
# https://docs.djangoproject.com/en/4.0/ref/settings/#default-auto-field

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'


# To view output in iframes
X_FRAME_OPTIONS = 'ALLOWALL'
XS_SHARING_ALLOWED_METHODS = ['POST','GET','OPTIONS', 'PUT', 'DELETE']

# LOGGING
logs_path = 'logs/debug.log'
os.makedirs(os.path.dirname(logs_path), exist_ok=True)
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'default': {
            'format': '[DJANGO] %(levelname)s %(asctime)s %(module)s '
                      '%(name)s.%(funcName)s:%(lineno)s: %(message)s'
        },
    },
    'handlers': {
        'file': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'formatter': 'default',
            'filename': logs_path,
        },

        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'default',
        }
    },
    'loggers': {
        '': {
            'handlers': ['console', 'file'],
            # 'handlers': ['file'],
            'level': 'DEBUG',
            'propagate': False,
        }
    },
}
# import logging
# logging.disable(logging.CRITICAL)


# Google OAuth Authentication
AUTHENTICATION_BACKENDS = [
    'django.contrib.auth.backends.ModelBackend',
    'allauth.account.auth_backends.AuthenticationBackend',
]

SOCIALACCOUNT_PROVIDERS = {
    'google': {
        'SCOPE': [
            'profile',
            'email',
        ],
        'AUTH_PARAMS': {
            'access_type': 'online',
        }
    }
}

# Login site redirect
if os.environ.get('DJANGO_SITE_ID') is None:
    SITE_ID = 0
else:
    # SITE_ID = 2 is what works for SPIN NERSC.
    SITE_ID = int(os.environ.get('DJANGO_SITE_ID'))    

LOGIN_REDIRECT_URL = '/retroapp'
LOGOUT_REDIRECT_URL = '/retroapp'