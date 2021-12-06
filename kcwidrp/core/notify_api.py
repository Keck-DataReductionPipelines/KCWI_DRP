import requests


def rti_config_data(config):

    url = config.rti_url
    user = config.rti_user
    pw = config.rti_pass,

    if type(pw) == tuple:
        try:
            pw = pw[0]
        except IndexError:
            pw = ''

    return url, user, pw


def send_data_api(data, config, logger):

    url, user, pw = rti_config_data(config)

    try:
        res = requests.get(url, params=data, auth=(user, pw))
        logger.info(f"Sending {data} to {res.request.url}")
    except requests.exceptions.RequestException as e:
        logger.error(f"Error caught while posting to {url}:")
        logger.error(e)
        return None

    return res

