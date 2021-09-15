#!/usr/bin/env python

import itertools
import axi

def print_board(board):
    print "\n".join(board)

def transpose(board):
    return [''.join(x) for x in zip(*board)]

def chunk(it, n):
    return zip(*[iter(it)]*n)

# 2 -> 00 01 10 11
def valid_rows(size=7, min_word=3):
    assert size < 30
    for i in xrange(1 << size):
        r = bin(i)[2:].zfill(size)
        if not all([x == '' or len(x) >= min_word for x in r.split("1")]):
            continue
        if not "0" in r:
            continue
        yield r.replace("1", "#").replace("0", ".")

def valid_boards(size=7, min_word=3):
    assert size % 2 == 1
    assert size < 30
    rows = lambda: valid_rows(size=size, min_word=min_word)
    for center in rows():
        if center != center[::-1]:
            continue
        for rs in itertools.product(rows(), repeat=size/2):
            board = list(rs) + [center] + [r[::-1] for r in rs][::-1]
            tboard = transpose(board)
            if "#" * size in tboard:
                continue
            if not all([x == '' or len(x) >= min_word for x in "#".join(tboard).split("#")]):
                continue
            yield tboard


def plot_board(board, x, y):
    paths = []
    #ratio = 0.4
    ratio = 1.0
    cell = ratio / len(board)
    #for dx in range(len(board) + 1):
    #    dx *= cell
    #    paths.append([(x + dx, y), (x + dx, y + ratio)])
    #for dy in range(len(board) + 1):
    #    dy *= cell
    #    paths.append([(x, y + dy), (x + ratio, y + dy)])

    for dx, row in enumerate(board):
        dx *= cell
        for dy, val in enumerate(row):
            dy *= cell
            if val == '.':
                continue
            paths.append([(x + dx, y + dy), (x + dx + cell, y + dy + cell)])
            paths.append([(x + dx + cell, y + dy), (x + dx, y + dy + cell)])

    return paths


def main():
    boards = list(valid_boards(size=7, min_word=3))
    if False:
        for i, vr in enumerate(boards):
            print i
            print_board(vr)
        print len(boards)

    assert len(boards) == 15 * 21 + 1

    if False:
        for y, brow in enumerate(chunk(boards[1:], 15)):
            for row in zip(*brow):
                print " ".join(row)
            print ""

    PADDING = 9 / 7.
    paths = []
    if True:
        for row in range(21):
            for dy in range(8):
                y = (row * PADDING + dy / 7.)
                paths.append([(0, y), (14 * PADDING + 1, y)])
        for col in range(15):
            for dx in range(8):
                x = (col * PADDING + dx / 7.)
                paths.append([(x, 0), (x, 20 * PADDING + 1) ])
        for y, brow in enumerate(chunk(boards[1:], 15)):
            for x, board in enumerate(brow):
                ps = plot_board(board, x * PADDING, y * PADDING)
                paths.extend(ps)

    d = axi.Drawing(paths)
    cli = axi.cli()
    cli.draw(d)

if __name__ == "__main__":
    main()
