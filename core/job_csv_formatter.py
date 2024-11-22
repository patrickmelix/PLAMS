import re

from scm.plams.core.basejob import Job
from scm.plams.core.enums import JobStatus
from scm.plams.core.logging import CSVFormatter


class JobCSVFormatter(CSVFormatter):

    def __init__(
        self,
        fmt=None,
        datefmt=None,
        style="%",
        validate=True,
        log_level=True,
        log_logger_name=False,
        log_time=True,
        log_lineno=False,
        counter_len=3,
        **csv_args,
    ):
        super().__init__(fmt, datefmt, style, validate, log_level, log_logger_name, log_time, log_lineno)
        self.counter_len = counter_len

    def format(self, record):
        if isinstance(record.msg, Job):
            record.msg = self.format_job(record.msg)
        return super().format(record)

    def format_job(self, job):
        message = {
            "job_base_name": job.name,
            "job_status": job.status,
            "job_path": "",
            "job_name": "",
            "job_ok": "",
            "job_check": "",
            "job_get_errormsg": "",
            "job_parent_name": "",
            "job_parent_path": "",
        }
        status_beginners = [JobStatus.CREATED, JobStatus.STARTED]
        status_intermediate = status_beginners.copy() + [JobStatus.REGISTERED, JobStatus.RUNNING]
        # status_final = (
        #     status_beginners.copy()
        #     + status_intermediate.copy()
        #     + [
        #         JobStatus.FINISHED,
        #         JobStatus.CRASHED,
        #         JobStatus.FAILED,
        #         JobStatus.SUCCESSFUL,
        #         JobStatus.COPIED,
        #         JobStatus.PREVIEW,
        #         JobStatus.DELETED,
        #     ]
        # )
        if job.status not in status_beginners:
            message.update({"job_path": job.path})
            message.update({"job_name": job.name})
            pattern = r"(\.\d{%i})+$" % self.counter_len
            base_name = re.sub(pattern, "", job.name)
            message.update({"job_base_name": base_name})
        if job.status not in status_intermediate:
            message.update({"job_ok": job.ok()})
            try:
                message.update({"job_check": job.check()})
            except TypeError:
                pass
            try:
                # this one it is not supported by the Job class but on many jobs they have it implemented!
                message.update({"job_get_errormsg": job.get_errormsg()})
            except (AttributeError, TypeError):
                pass
        if job.parent:
            add = {
                "job_parent_name": job.parent.name,
                "job_parent_path": job.parent.path,
            }
            message.update(add)
        return message
