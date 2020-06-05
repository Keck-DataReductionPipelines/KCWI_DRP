from keckdrpframework.config.framework_config import ConfigClass
import logging.config
import pkg_resources

pkg = 'kcwidrp'


def test_parse_kcwi_config():

    kcwi_config_file = 'configs/kcwi.cfg'
    kcwi_config_fullpath = pkg_resources.resource_filename(
        pkg, kcwi_config_file)
    kcwi_config = ConfigClass(kcwi_config_fullpath, default_section='KCWI')


def test_parse_log_config():

    framework_logcfg_file = 'configs/logger.cfg'
    framework_logcfg_fullpath = \
        pkg_resources.resource_filename(pkg, framework_logcfg_file)
    # we are only checking that the config file can be parsed,
    # not that the log file can be created
    try:
        logging.config.fileConfig(framework_logcfg_fullpath)
    except FileNotFoundError:
        pass


def test_parse_framework_config():

    framework_config_file = 'configs/framework.cfg'
    framework_config_fullpath = pkg_resources.resource_filename(
        pkg, framework_config_file)
    framework_config = ConfigClass(framework_config_fullpath)