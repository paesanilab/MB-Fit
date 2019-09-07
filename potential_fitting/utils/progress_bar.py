from . import system


class ProgressBar(object):

    def __init__(self, start=0, end=100):
        self.start = start
        self.end = end

    def update(self, progress):
        percent = int(100 * (progress - self.start) / (self.end - self.start))
        string = "|" + "="*percent + " "*(100 - percent) + "|"

        system.format_print(string, bold=True, color=system.Color.BLUE, replace_with_next_line=True)

    def finish(self):

        percent = 100
        string = "|" + "="*percent + " "*(100 - percent) + "|"

        system.format_print(string, bold=True, color=system.Color.GREEN, replace_with_next_line=True)

        print()
