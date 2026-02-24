

from pathlib import Path


def pytest_addoption(parser):
    parser.addoption(
        '--reference',action='store',default=None,
        help='Folder with reference solution.'
        )


def pytest_sessionstart(session):
    reference = session.config.getoption('--reference')

    if reference is not None:
        if not Path(reference).is_dir():
            raise ValueError(f'Reference folder {reference} does not exist!')
        session.config.reference = reference
        print('- Selected reference folder:',reference)
    