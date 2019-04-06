import itertools

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
    path = ""
    ratio = 0.4
    cell = ratio / len(board)
    x *= ratio * 1.1
    y *= ratio * 1.1
    for dx in range(len(board) + 1):
        dx *= cell
        path += "M {} {} v {} ".format(x + dx, y, ratio)
    for dy in range(len(board) + 1):
        dy *= cell
        path += "M {} {} h {} ".format(x, y + dy, ratio)

    for dx, row in enumerate(board):
        dx *= cell
        for dy, val in enumerate(row):
            dy *= cell
            if val == '.':
                continue
            path += "M {} {} l {} {} ".format(x + dx, y + dy, cell, cell)
            path += "m -{} 0 l {} -{} ".format(cell, cell, cell)


    style = "fill:none; stroke:#000000; stroke-width:0.01; stroke-opacity:1"
    print '<g style="{}"><path d="{}" /></g>'.format(style, path)

def svg_start():
    print """
<svg
    xmlns="http://www.w3.org/2000/svg"
    width="8.5in"
    height="11in"
    viewBox="0 0 8.5 11">
<g id="layer1">
"""

def svg_end():
    print "</g></svg>"


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

    if True:
        svg_start()
        for y, brow in enumerate(chunk(boards[1:], 15)):
            for x, board in enumerate(brow):
                plot_board(board, x, y)
        svg_end()


if __name__ == "__main__":
    main()
