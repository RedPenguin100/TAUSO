import os
import pathlib


class FileUtil(object):
    root_dir = None

    @classmethod
    def set_root_dir(cls):
        script_directory = pathlib.Path(__file__).parent.resolve()
        # root = str(script_directory.parents[1])
        root = script_directory
        cls.root_dir = root

    @classmethod
    def get_root_dir(cls):
        if FileUtil.root_dir is None:
            cls.set_root_dir()
        return cls.root_dir

    @classmethod
    def get_output_dir(cls):
        return '/dev/shm/tauso/_raccess/output'

    @classmethod
    def get_output_path(cls, file_name):
        file_dir = cls.get_output_dir()
        file_path = os.path.join(file_dir, file_name)
        return file_path
