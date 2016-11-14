"""
This module aggregates all the required information associated with a job list file
Plays the role of "model" in MVC
"""

import subprocess, re, os
from datetime import datetime, timedelta


###################################################################################################
# CLASSES
###################################################################################################

# Data container for the information associated with a job
class JobLog:
    # These are the fields in a JobLog
    __slots__ = ['job_id', 'job_full_id', 'job_name', 'job_dependencies',
                 'status', 'exit_status', 'MUGQIC_exit_status',
                 'walltime', 'start_date', 'end_date', 'cput', 'cpu_to_real_time_ratio',
                 'mem', 'vmem', 'vmem_to_mem_ratio', 'limits', 'queue',
                 'username', 'group', 'nodes', 'path']
    # You can initialize a JobLog via keyword arguments, eg:
    # JobLog(job_id='123', limits='ppn=2', mem='40Gb')
    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    # Useful for debugging purposese
    def __repr__(self):
        try:
            return self.job_name + ' (' + self.job_id + ')'
        except:
            return 'Undefined JobLog instance'


class JobDependencies:
    """
    Data container for the dependencies of a job

    You can iterate over this container, or call 'str' on it to get a colon-delimited string
    """
    def __init__(self, dependency_string):
        """
        :param dependency_string: str from job list file, eg '123:456'
        """
        # Eg. '123:456' -> ['123', '456']
        self.dependencies = [id for id in dependency_string.split(':') if id != '']

    def __len__(self):
        return len(self.dependencies)

    # This allows us to iterate over instances of this class
    # Eg. "for id in dependencies:..."
    def __iter__(self):
        return iter(self.dependencies)

    # Eg: '123:456:987'
    def __str__(self):
        return ':'.join(str(dep) for dep in self.dependencies)


class RE:
    """
    Wrapper for perl-like regex matching

    Eg:
    if RE.match(<pattern>, <text>):
        do_something(RE.group(1), RE.group(1))
    """
    @classmethod
    def search(cls, pattern, string, flag=0):
        m = re.search(pattern, string, flag)
        if m != None:
            cls.groups = m.groups()
            return True

        return False

    @classmethod
    def group(cls, num):
        """
        :return: the 'num'th expression captured in parentheses
        """
        return cls.groups[num - 1]


class MemorySize:
    """
    Class to keep track of memory sizes (eg. '123kb') and be useful for displaying

    Use the 'bytes' int field
    """
    def __init__(self, memory_size):
        """
        If no unit is given, bytes are assumed

        :param memory_size: formatted string (eg. '2M')
        """
        if RE.search('^(\d*\.?\d*)g', memory_size, re.IGNORECASE):
            self.bytes = int(float(RE.group(1)) * 1024 ** 3)
        elif RE.search('^(\d*\.?\d*)m', memory_size, re.IGNORECASE):
            self.bytes = int(float(RE.group(1)) * 1024 ** 2)
        elif RE.search('^(\d*\.?\d*)k', memory_size, re.IGNORECASE):
            self.bytes = int(float(RE.group(1)) * 1024)
        elif RE.search('^(\d+)b', memory_size, re.IGNORECASE):
            self.bytes = int(float(RE.group(1)))
        else:
            raise Exception('Can\'t convert ' + memory_size + ' to a value in bytes')

    # The default representation of a memory size in Gb
    def __str__(self):
        return '{0:.2f}'.format(float(self.bytes) / 1024 ** 3) + ' GiB'

    # We must define the '<' so we can sort the MemorySizes later
    def __lt__(self, other):
        return self.bytes < other.bytes


###################################################################################################
# HELPERS
###################################################################################################

def conditional_assign(obj, attribute_name, value):
    """
    Assign obj.attribute_name = value if the attribute is not defined
    """
    if not hasattr(obj, attribute_name):
        setattr(obj, attribute_name, value)


def run_command(cmd_list):
    """
    E.g. ['showjobs', '-o', '-u', 'nreinhardt']
    :param cmd_list: list of command, options, and values
    :return: stdout of command results, with leading and trailing whitespace removed
    :raises: Exception if command fails
    """
    # DEBUG
    # print('Running command: ' + ' '.join(cmd_list))

    process = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()
    ret_code = process.wait()

    if ret_code != 0 or output.strip() == '':
        raise Exception('Command ' + str(cmd_list) + ' failed or did not return output')

    # If the output isn't a string, decode it
    if not isinstance(output, str):
        output = output.decode('utf-8')

    return output.strip()


def get_file_date(filename):
    """
    :return: seconds since epoch of date file status last changed, or None if the date can't be found
    """
    try:
        return int(os.path.getctime(filename))
    except:
        return None


def get_user():
    try:
        return run_command(['whoami'])
    except:
        return None


def percent(num):
    """
    Format num as a percent
    :return: str
    """
    return '{0:.0f}%'.format(100 * num)


###################################################################################################
# JOB LIST FILE
###################################################################################################

def parse_job_list_file(filename):
    """
    :param filename: path to job list file
    :return: list of JobLogs
    """
    job_logs = []

    with open(filename) as f:
        for line in f.readlines():
            id, name, dependency_string, log_file_path = line.strip().split('\t')

            dependencies = JobDependencies(dependency_string)

            # Get the directory of the job list file, and then add on the relative path in the job list file
            # to get the full path to the log file
            full_log_file_path = os.path.abspath(os.path.join(os.path.dirname(filename), log_file_path))

            job_logs.append(JobLog(job_name=name, job_id=id, job_dependencies=dependencies, path=full_log_file_path))

    return job_logs


###################################################################################################
# CLUSTER JOB LOG
###################################################################################################

def parse_cluster_job_log_paths(job_output_dir, job_logs):
    """
    Get the MUGQICexitStatus from each job's full log

    :param job_output_dir: path to directory containing all the job logs
    :param job_logs: list of JobLogs
    :return: modified list of JobLogs
    """
    for log in job_logs:
        if os.path.exists(log.path):
            # If the output file exists, then the job must have ended at the date the file was created
            log.end_date = datetime.fromtimestamp(get_file_date(log.path))

            with open(log.path) as f:
                for line in f.readlines():
                    if RE.search('MUGQICexitStatus:(\S+)', line):
                        status = RE.group(1)

                        # Assign status based on the MUGQIC exit status
                        log.MUGQIC_exit_status = status
                        log.status = 'SUCCESS' if status == '0' else 'FAILED'

    return job_logs


###################################################################################################
# FILTER BY SUCCESS
###################################################################################################

def filter_by_success(job_logs, keep_successful, keep_unsuccessful):
    """
    :param job_logs: list of JobLogs
    :param keep_successful: whether to keep successful jobs
    :param keep_unsuccessful: whether to keep jobs that have not succeeded yet (blocked or failed)
    :return: only those JobLogs that search the given criteria
    """
    is_successful = lambda log: hasattr(log, 'MUGQIC_exit_status')\
                                and log.MUGQIC_exit_status == '0'

    if keep_successful:
        job_logs = [log for log in job_logs if is_successful(log)]
    elif keep_unsuccessful:
        job_logs = [log for log in job_logs if not is_successful(log)]

    return job_logs


def assign_cpu_to_real_time_ratio(job_logs):
    """
    Assign "cput / walltime" ratio to each JobLog if possible
    :param job_logs: list of JobLogs
    :return: modified list of JobLogs
    """
    for log in job_logs:
        if hasattr(log, 'cput') and hasattr(log, 'walltime') and log.walltime.total_seconds() != 0:
            log.cpu_to_real_time_ratio = percent(log.cput.total_seconds() / log.walltime.total_seconds())

    return job_logs


def assign_vmem_to_mem_ratio(job_logs):
    """
    Assign "vmem / mem" ratio to each JobLog if possible
    :param job_logs: list of JobLogs
    :return: modified list of JobLogs
    """
    for log in job_logs:
        if hasattr(log, 'vmem') and hasattr(log, 'mem') and log.mem.bytes != 0:
            log.vmem_to_mem_ratio = percent(float(log.vmem.bytes) / float(log.mem.bytes) - 1)

    return job_logs


###################################################################################################
# SHOWJOBS PARSING
###################################################################################################

def get_end_date(job_log_list):
    """
    :param job_log_list: list of JobLogs
    :return: the maximum end date of all jobs in format YYYY-MM-DD, or None if we can't find it
    """
    try:
        # The end dates of all the jobs in Unix timestamps
        # Filter out all 'None' values when we can't find a date
        end_date_timestamps = filter(None, [get_file_date(log.path) for log in job_log_list])

        # Return the maximum date
        end_date = max(datetime.fromtimestamp(stamp) for stamp in end_date_timestamps)
        return end_date.strftime('%Y-%m-%d')
    except:
        return None


def get_all_showjobs_output(job_list_filename, job_log_list):
    """
    :param job_list_filename: path to the job list file
    :param job_log_list: list of JobLogs
    :return: string containing all possible showjobs output for this job log
    """
    # Get the end date of the latest job
    query_start_date = datetime.fromtimestamp(int(get_file_date(job_list_filename))).strftime('%Y-%m-%d')
    # Options must be passed as a list of strings
    start_date_option = ['-s', query_start_date] if query_start_date else []

    # Get the earliest start date
    query_end_date = get_end_date(job_log_list)
    end_date_option = ['-e', query_end_date] if query_end_date else []

    user = get_user()
    user_option = ['-u', str(user)] if user else []

    # Run a 'showjobs' query on the specified date range for this user
    try:
        showjobs_results = run_command(['showjobs'] + start_date_option + end_date_option + user_option)
    except:
        return None

    return showjobs_results


def get_showjobs_output(job_list_filename, job_logs):
    """
    :param job_list_filename: path to the job list file
    :param job_logs: list of JobLogs
    :return: dict mapping from job id -> list of lines of relevant showjobs results
    """
    id_to_entry = dict()

    all_results = get_all_showjobs_output(job_list_filename, job_logs)
    if all_results == None:
        return None

    ids = set([log.job_id for log in job_logs])

    for entry in all_results.split('-' * 80):
        entry_lines = entry.splitlines()

        for line in entry_lines:
            if RE.search('Job Id *: (\d+)', line): # We found the job id line
                id = RE.group(1)

                # Don't save this output entry unless this id is one of the ids we are looking for
                if id in ids:
                    id_to_entry[id] = entry_lines
                break

    return id_to_entry


def parse_showjobs_output(id_to_showjobs_entry, job_log):
    """
    Add fields to job_log by parsing the showjobs output

    :param id_to_showjobs_entry: dict mapping from job id -> list of showjobs output lines
    :return: modified JobLog
    """
    # Make sure that we have showjobs output for this job
    if job_log.job_id in id_to_showjobs_entry.keys():
        # Get the relevant part of the showjobs output for this job
        entry = id_to_showjobs_entry[job_log.job_id]

        for line in entry:
            if RE.search('^Job Id *: (\d+)$', line):
                job_log.job_full_id = RE.group(1)
            elif RE.search('^Start Time *: (.*)$', line):
                job_log.start_date = datetime.strptime(RE.group(1), '%a %b %d %X %Y')
            elif RE.search('^User Name *: (\S+)$', line):
                job_log.username = RE.group(1)
            elif RE.search('^Group Name *: (\S+)$', line):
                job_log.group = RE.group(1)
            elif RE.search('^CPUTime *: (\d+):(\d+):(\d+)$', line):
                job_log.cput = timedelta(hours=int(RE.group(1)), minutes=int(RE.group(2)), seconds=int(RE.group(3)))
            elif RE.search('^Memory Used *: (\S+)$', line):
                job_log.mem = MemorySize(RE.group(1))
            elif RE.search('^vmem Used *: (\S+)$', line):
                job_log.vmem = MemorySize(RE.group(1))
            elif RE.search('^Wallclock Duration *: (\d+):(\d+):(\d+)$', line):
                job_log.walltime = timedelta(hours=int(RE.group(1)), minutes=int(RE.group(2)), seconds=int(RE.group(3)))
            elif RE.search('^Queue Name *: (\S+)$', line):
                job_log.queue = RE.group(1)
            elif RE.search('^Exit Code *: (\S+)$', line):
                status = RE.group(1)
                job_log.exit_status = status
                # status is SUCCESS if and only if exit_status == '0'
                if status == '0':
                    conditional_assign(job_log, 'status', 'SUCCESS')
                else:
                    conditional_assign(job_log, 'status', 'FAILED')
            elif RE.search('^End Time *: (.*)$', line):
                job_log.end_date = datetime.strptime(RE.group(1), '%a %b %d %X %Y')

    return job_log


###################################################################################################
# CHECKJOB PARSING
###################################################################################################

def get_checkjob_output(job_id):
    '''
    :return: list of lines containing checkjob output for this job
    '''
    try:
        return run_command(['checkjob', '-v', str(job_id)]).splitlines()
    except:
        return ''


def parse_checkjob_output(checkjob_results, job_log):
    """
    Add fields to job_log by parsing the showjobs output
    Only assign fields if they have not been defined already (don't overwrite)

    :param checkjob_results: list of lines containing checkjob output for this job
    :return: modified job_log
    """
    for index, line in enumerate(checkjob_results):
        if RE.search('^job (\d+)$', line):
            conditional_assign(job_log, 'job_full_id', RE.group(1))
        elif RE.search('^State: (\S+)', line):
            state = RE.group(1)
            if state == 'Idle':
                conditional_assign(job_log, 'status', 'INACTIVE')
            elif state == 'Running':
                conditional_assign(job_log, 'status', 'ACTIVE')
            elif state == 'Removed':
                conditional_assign(job_log, 'status', 'FAILED')
        elif RE.search('^Completion Code: (\d+) +Time: (.*)$', line):
            conditional_assign(job_log, 'exit_status', RE.group(1))
            # We have to add in the current year manually
            conditional_assign(job_log, 'end_date', datetime.strptime(RE.group(2) + ' ' + str(datetime.now().year), '%a %b %d %H:%M:%S %Y'))
        elif RE.search('^Creds: *user:(\S+) *group:(\S+) *class:(\S+)', line):
            conditional_assign(job_log, 'username', RE.group(1))
            conditional_assign(job_log, 'group', RE.group(2))
        # Only assign the walltime if the job has finished
        elif RE.search('^WallTime:\s+(\d+):(\d+):(\d+)', line) and hasattr(job_log, 'MUGQIC_exit_status'):
            conditional_assign(job_log, 'walltime', timedelta(hours=int(RE.group(1)), minutes=int(RE.group(2)), seconds=int(RE.group(3))))
        elif RE.search('^StartTime: (.*)$', line):
            conditional_assign(job_log, 'start_date', datetime.strptime(RE.group(1) + ' ' + str(datetime.now().year), '%a %b %d %H:%M:%S %Y'))
        elif RE.search('^Allocated Nodes:', line):
            # This field's value is on the following line, for example
            # Allocated Nodes:
            # [r2a - 1:3][r3a - 1:3][r4a - 1:3]
            next_line = checkjob_results[index + 1]
            conditional_assign(job_log, 'nodes', next_line)
        elif RE.search('^OuptutFile: *(\S+):(//\S+)$', line):
            conditional_assign(job_log, 'path', os.path.abspath(RE.group(1)))
        elif RE.search('^Submit Args:', line):
            # Try to match each limit pattern
            for pattern in ['(nodes:\d+)', '(ppn=\d+)', '(walltime=\d+:\d+:\d+)',
                            '(vmem=\d+\S)', '(mem=\d+\S)']:
                # If we have found a resource limitation, add to a colon-delimited text
                # Eg. ppn=1:walltime=24:00:00:vmem=10g
                if RE.search(pattern, line):
                    if hasattr(job_log, 'limits'):
                        job_log.limits += ':' + RE.group(1)
                    else:
                        job_log.limits = RE.group(1)

    return job_log


###################################################################################################
# CREATE JOB LOGS
###################################################################################################

def create_job_logs(job_list_file, keep_successful, keep_unsuccessful, minimal_detail=False):
    """
    Aggregate information about each job and return a list of JobLogs

    If the information for a field cannot be found, then the field is not defined
    :param job_list_file: path to the job list file
    :param keep_successful: whether to keep the successful jobs
    :param keep_unsuccessful:  whether to keep the unsuccessful or blocked jobs
    :return: list of JobLogs
    """
    job_logs = parse_job_list_file(job_list_file)

    # Path to the directory containing all the job logs
    job_output_dir = os.path.dirname(job_list_file)

    # Go into each job's full log, and extract the results
    job_logs = parse_cluster_job_log_paths(job_output_dir, job_logs)

    if not minimal_detail:
        # Remove the successful/unsuccessful jobs as dictated by the given options
        job_logs = filter_by_success(job_logs, keep_successful, keep_unsuccessful)

        # dict mapping from id -> list of lines of relevant showjobs results
        id_to_showjobs_entry = get_showjobs_output(job_list_file, job_logs)

        for log in job_logs:
            # Fill in fields using information from showjobs
            if id_to_showjobs_entry:
                log = parse_showjobs_output(id_to_showjobs_entry, log)

            # Fill in fields using information from checkjob
            checkjob_results = get_checkjob_output(log.job_id)
            log = parse_checkjob_output(checkjob_results, log)

        # TODO: more fields can be gathered from the 'qstat' command (eg. session, account)

        job_logs = assign_cpu_to_real_time_ratio(job_logs)
        job_logs = assign_vmem_to_mem_ratio(job_logs)


    return job_logs