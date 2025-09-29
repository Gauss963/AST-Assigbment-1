from typing import TextIO, Union


def headline(id_char: Union[str, int], iread: TextIO) -> None:
    if isinstance(id_char, int):
        if not (0 <= id_char <= 0x10FFFF):
            raise ValueError("id_char (int) out of valid Unicode range.")
        target = chr(id_char)
    elif isinstance(id_char, str) and len(id_char) == 1:
        target = id_char
    else:
        raise ValueError("id_char must be a single-character string or an integer code point.")

    while True:
        line = iread.readline()
        if line == "":
            return
        if line and line[0] == target:
            return