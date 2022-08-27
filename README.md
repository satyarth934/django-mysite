# django-mysite
![example workflow](https://github.com/satyarth934/django-mysite/actions/workflows/docker-image.yml/badge.svg)


# Rancher/Spin launch instructions
1. Finalize the Django code.
2. Push to GitHub
    1. Only for safe keeping.
    2. Later we can explore the upload to [registry.nersc.gov](http://registry.nersc.gov) in GitHub actions.
3. Docker build and push to [registry.nersc.gov](http://registry.nersc.gov)
    1. It is recommended to do this on a linux machine.
5. Update the website URL on the Spin Workload instance.
6. Open Shell from Spin for the website instance.
    1. Navigate to the django project root dir (should be there by default as it's specified in the Dockerfile).
    2. Run the following command to launch the server:
    `python3 manage.py runserver 0.0.0.0:8000`
6. Open this link in a separate browser tab to open the webpages and webapps:
[http://molinv.molinv.development.svc.spin.nersc.org/molinv/](http://molinv.molinv.development.svc.spin.nersc.org/molinv/)
