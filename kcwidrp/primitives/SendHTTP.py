import os
import requests
import time

from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
                                                    kcwi_fits_reader


class SendHTTP(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
    
    def _perform(self):

        if not self.action.args.koaid:
            self.logger.error(f"Encountered a file with no KOA ID: {self.action.args.name}")
            return self.action.args
        
        self.logger.info(f"Alerting RTI that {self.action.args.name} is ready for ingestion")

        url = self.config.instrument.rtiurl
        data = {
            'instrument': 'KCWI',
            'koaid': self.action.args.koaid,
            'ingesttype': 'lev2',
            'datadir': str(self.config.instrument.output_directory),
            'start': str(self.action.args.ingest_time),
            'reingest': True,
            'testonly': True,
            'dev': True
        }
        
        attempts = 0
        limit = self.config.instrument.rti_attempts
        while attempts < limit:
            post = self.post_url(url, data)
            if post is None:
                t = self.config.instrument.rti_retry_time
                attempts += 1
                self.logger.error(f"Waiting {t} seconds to attempt again... ({attempts}/{limit})")
                time.sleep(t)
            else:
                break


        self.logger.info(f"Post returned status code {post.status_code}")

        return self.action.args
    
    def post_url(self, url, data):
        try:
            self.logger.info(f"Posting to RTI with KOAID {data['koaid']}")
            post = requests.post(url, data = data)
        except requests.exceptions.RequestException as e:
            self.logger.error(f"Error caught while posting to {url}:")
            self.logger.error(e)
            return None
        return post
