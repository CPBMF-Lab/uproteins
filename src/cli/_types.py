# Copyright © 2025 Eduardo Vieira de Souza
# Copyright © 2025 Adriana Canedo
# Copyright © 2025 Cristiano Valim Bizarro
# Copyright © 2025 Bruno Maestri A Becker
#
# This file is part of uProteInS.
#
# uProteInS is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# uProteInS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# uProteInS. If not, see <https://www.gnu.org/licenses/>.


"""Utilitarian module containing a number of helper functions designed to serve
as values for the :param:`type` parameter of the :mod:`argparse` module.

Those functions raise :exc:`ArgumentError` when the str value passed to them
fails to conform to their expected format, else convert the value to the
pertinent type and return it. Pathes are returned in their absolute forms.

All functions can receive a None value, in which case they just return None.
"""

import argparse
import pathlib
import shutil
import re
import typing as t


T = t.TypeVar('T')


class Executable:
    def __init__(self, default: str) -> None:
        """Executable type. This type is a callable class. Initialize it at the
        type field with the default value. The type checking won't fail even
        if it's not an Executable.

        When called
        -----------
        Receive a str path and raise an :exc:`ArgumentError` if the value
        is not a valid executable, else return a :obj:`Path` object.

        Create an :obj:`Executable` with `'exec'` as default.
        Passing 'exec' won't raise even if `'exec'` is not an executable.
        >>> exec = Executable('exec')
        >>> exec('exec')
        Path('path/to/exec')

        But passing another invalid executable raises.
        >>> exec('foo')
        ArgumentTypeError: invalid Executable value: 'foo'

        If the arg is a valid executable, it doesn't raise.
        >>> exec('valid')
        Path('path/to/valid.exe')
        """
        self.default = default

    def __call__(self, val: str) -> pathlib.Path:
        """Receive a str path and raise an :exc:`ArgumentError` if the value
        is not a valid executable, else return a :obj:`Path` object.
        """
        if val == self.default:
            return pathlib.Path(val)

        which = shutil.which(val)
        if which is None:
            # Necessary to raise ArgumentTypeError with custom message for
            # correct display by argparse
            raise argparse.ArgumentTypeError(
                f"invalid Executable value: '{val}'"
            )

        return pathlib.Path(which).absolute()


def FilePath(val: str) -> pathlib.Path:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    the path to a file, else return a :obj:`Path` object.
    """
    path = pathlib.Path(val)
    if not (path.exists() and path.is_file()):
        raise TypeError
    # Test if it's readable
    with open(path, 'r'):
        pass
    return path.absolute()


def FileName(val: str) -> pathlib.Path:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a filename, else return a :obj:`Path` object.

    A filename is considered to be a path to a file or a path that doesn't yet
    exist.
    """
    path = pathlib.Path(val)
    if path.exists() and not path.is_file():
        raise TypeError
    if path.exists():
        # Test if it's writable
        with open(path, 'a'):
            pass
    return path.absolute()


def DirectoryPath(val: str) -> pathlib.Path:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a valid directory, else return a :obj:`Path` object.
    """
    path = pathlib.Path(val)
    if not (path.exists() and path.is_dir()):
        raise TypeError
    return path.absolute()


def DirectoryName(val: str) -> pathlib.Path:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a valid directoryname, else return a :obj:`Path` object.

    A valid directoryname is considered to be a path to a directory or a path
    that doesn't yet exist.
    """
    path = pathlib.Path(val)
    if path.exists() and not path.is_dir():
        raise TypeError
    return path.absolute()


def Codon(val: str) -> str:
    """Receive a str val and raise an :exc:`ArgumentError` if the value is not
    a valid codon, else return the str untouched.

    A valid codon is considered to be any three letter str composed only of
    the `ATCG` characteres.
    """
    if val is None:
        return val

    if re.fullmatch(r'[ATCG]{3}', val) is None:
        raise TypeError
    return val


def PositiveInt(val: str) -> int:
    n = int(val)
    if n < 1:
        raise TypeError
    return n


def Memory(val: str) -> str:
    if re.fullmatch(r'[0-9]*[1-9]+[0-9]*[kKgGmM]', val) is None:
        raise TypeError
    return val


class YesOrNoBooleanAction(argparse.Action):
    def __init__(
        self,
        option_strings,
        dest,
        default=False,
        required=False,
        help=None,
        metavar=None,
    ):
        if not isinstance(default, bool):
            raise TypeError(f'default parameter must be a bool, not {default}')

        if metavar is None:
            # Set metavar to capitalized default option
            metavar = '{Y/n}' if default else '{y/N}'

        super().__init__(
            option_strings,
            dest,
            default=default,
            required=required,
            help=help,
            metavar=metavar,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        if 'YES'.startswith(values.upper()):   # pyright: ignore
            setattr(namespace, self.dest, True)
        elif 'NO'.startswith(values.upper()):  # pyright: ignore
            setattr(namespace, self.dest, False)
        else:
            parser.error(
                f"argument {self.option_strings[0]}: invalid choice: "
                f"'{values}' (choose from 'YES','NO')"
            )


class CommaListAction(argparse.Action):
    """Comma list action. This will generate a Python list and apply
    the type used in initialization to each value.

    When used
    ---------
    Receive a str containing a comma-separated list of values. For each value,
    apply the type conversion with proper error handling. Return a Python list.

    >>> parser.add_argument('nums', type=int, action=CommaListAction)
    >>> args = parser.parse('1,2,3')
    >>> args.nums
    [1, 2, 3]
    >>> args = parser.parse('1,2,dd')
    ArgumentTypeError: invalid int value: 'dd'
    """
    def __init__(
        self,
        option_strings,
        dest,
        nargs=None,
        const=None,
        default=None,
        type=None,
        required=False,
        help=None,
        metavar=None,
    ):
        if type is None:
            type = lambda x: x

        super().__init__(
            option_strings,
            dest,
            nargs,
            const,
            default=default,
            type=self._Applicative(type),
            required=required,
            help=help,
            metavar=metavar,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)

    def format_usage(self):
        return f'{self.metavar},[{self.metavar}...]'

    class _Applicative(t.Generic[T]):
        def __init__(self, type: t.Callable[[str], T]):
            self.type = type

        def __call__(self, val: str) -> list[T]:
            return_list: list[T] = []
            for to_be_converted in val.split(','):
                try:
                    converted = self.type(to_be_converted)
                    return_list.append(converted)
                # This is necessary for the correct error display by argparse
                except (ValueError, TypeError):
                    raise argparse.ArgumentTypeError(
                        "invalid %s value: '%s'"
                        % (self.type.__name__, to_be_converted)
                    )

            return return_list
