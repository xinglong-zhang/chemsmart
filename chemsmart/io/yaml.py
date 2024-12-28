import logging
import os

import yaml

from chemsmart.utils.mixins import YAMLFileMixin

logger = logging.getLogger(__name__)


class YAMLFile(YAMLFileMixin):
    def __init__(self, filename):
        self.filename = filename
