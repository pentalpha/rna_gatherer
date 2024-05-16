import os
import sys
from gatherer.util import runCommand

class Pipeline:
    def __init__(self, args, confs, funcs, output_dir):
        self.nSteps = len(funcs)
        self.outputdir = os.path.abspath(output_dir)

        self.args = args
        self.confs = confs

        self.stepFuncs = {}
        self.order = list()
        self.stepDir = {}
        self.tempDir = {}
        i = 0
        for pair in funcs:
            self.stepFuncs[pair[0]] = pair[1]
            self.order.append(pair[0])
            dir_name = self.outputdir + "/step_" + str(i + 1) + "-" + pair[0]
            print(dir_name)
            self.stepDir[pair[0]] = dir_name
            self.tempDir[pair[0]] = dir_name + "-tmp"
            i += 1
        if not os.path.exists(self.outputdir):
            print("Creating " + self.outputdir)
            code = runCommand("mkdir " + self.outputdir)
            if code != 0:
                print("No permission to create " + self.outputdir + ", cannot procede.")
                self.ready = False
            else:
                self.ready = True
        else:
            self.ready = True

    def get_step_name(self, x):
        return self.order[x]

    def get_step_order(self, x):
        return self.order.index(x)

    def get_dir(self, x):
        return self.stepDir[self.get_step_name(x)]

    def run(self, start_from = "-1", stop_at = "-1"):
        print("Running pipeline")
        try:
            int(start_from)
            start_from = int(start_from)
        except ValueError:
            start_from = self.get_step_order(start_from) + 1

        try:
            int(stop_at)
            stop_at = int(stop_at)
        except ValueError:
            stop_at = self.get_step_order(stop_at)

        if not self.ready:
            print("Not ready to start pipeline.")
            return

        startingStep = 1
        if start_from > 0:
            startingStep = start_from
            if startingStep > self.nSteps:
                sys.exit("This step does not exist")
            elif startingStep > 1:
                if not os.path.exists(self.get_dir(startingStep-2)):
                    sys.exit("The previous step to Step " + str(startingStep) + " has not been done yet.")
                    print("Starting from " + str(startingStep))
        else:
            for name, path in self.stepDir.items():
                if os.path.exists(path):
                    print("Skipping step " + str(startingStep))
                    startingStep += 1

        #running necessary steps
        limit = len(self.stepFuncs)
        if stop_at > 0:
            stopAt = stop_at
            if stopAt > 0:
                limit = stopAt
                print("Stoping at " + str(limit))
        #print("entering steps loop from " + str(startingStep-1) + " " + str(limit))
        #print(str(range(startingStep-1, limit)))
        for i in range(startingStep-1, limit):
            print(str(i))
            step = self.get_step_name(i)
            #print("entered")
            print("--- STEP " + str(i+1) + ": " + step + " ---")

            #create temporary dir to store files from next step
            if os.path.exists(self.tempDir[step]):
                runCommand("rm -Rf " + self.tempDir[step])
            runCommand("mkdir " + self.tempDir[step])
            print(self.tempDir[step])
            #run step
            success = self.stepFuncs[step](self.args, self.confs, self.tempDir[step], self.stepDir)

            if success:
                #move results from temporary dir to permanent one
                if os.path.exists(self.stepDir[step]):
                    runCommand("rm -Rf " + self.stepDir[step])
                runCommand("mv " + self.tempDir[step] + " " + self.stepDir[step])
            else:
                print("Step " + str(i+1) + " was not successful.")
                break
