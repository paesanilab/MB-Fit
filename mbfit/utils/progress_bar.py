from . import system


class ProgressBar(object):

    def __init__(self, start=0, end=100):
        self.start = start
        self.end = end
        self.progression = self.start

    def progress(self, amount):
        self.progression += amount
        self.update(self.progression)

    def update(self, progress):
        self.progression = progress
        percent = int(100 * (self.progression - self.start) / (self.end - self.start))
        string = "|" + "="*percent + " "*(100 - percent) + "|"

        system.format_print(string, bold=True, color=system.Color.BLUE, replace_with_next_line=True)

    def finish(self):
        self.progression = self.end

        percent = 100
        string = "|" + "="*percent + " "*(100 - percent) + "|"

        system.format_print(string, bold=True, color=system.Color.GREEN, replace_with_next_line=True)

        print()
