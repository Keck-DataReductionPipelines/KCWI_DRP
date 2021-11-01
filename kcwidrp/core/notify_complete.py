import requests


def rti_config_data(config):

    url = config.rti_url
    user = config.rti_user
    pw = config.rti_pass,

    return url, user, pw


def send_api_complete(config, utdate, logger):

    data = {
        'instrument': 'KCWI',
        'utdate': utdate,
        'ingesttype': 'lev2'
    }

    url, user, pw = rti_config_data(config)

    try:
        res = requests.get(url, params=data, auth=(user, pw))
        logger.info(f"Sending Complete Status {res.request.url}")
    except requests.exceptions.RequestException as e:
        logger.error(f"Error caught while posting to {url}:")
        logger.error(e)
        return None

    return res

