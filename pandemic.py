""" 
Scarce resources are assigned as cases come in and cases come in in the order demographics are listed.
If a vulnerable group is more likely to be given the resources, they should be listed first.
"""
#%%

#%%
class Population:
    def __init__(
        self, r0, incubation, duration, asymptomatic, demographics, interventions, resources
    ):
        self.day = 0
        self.pop = 0
        self.infected = 0
        self.dead = 0
        self.r0 = r0
        self.inc = incubation
        self.dur = duration
        self.asymp = asymptomatic
        self.sick_dur = duration - incubation
        self.demographics = []
        self.intervention_demos = {}
        self.interventions = []
        self.resources = []
        self.inf_hist = []
        self.death_hist = []

        for d in demographics:
            self.pop += int(d["size"])
            self.infected += int(d["infected"])
            self.demographics.append(
                self.Demographic(
                    d["name"], int(d["size"]), int(d["infected"]), float(d["cfr"]), duration
                )
            )

        for i in interventions:
            if i["name"]:
                self.interventions.append(
                    self.Intervention(
                        i["name"],
                        i["criterion"],
                        i["threshold"],
                        i["demo_info"],
                        i["variolation"],
                        i["quarantine"],
                    )
                )

        for r in resources:
            if r["name"]:
                self.resources.append(self.Resource(r["name"], r["number"], r["demo_info"]))

    def init_variolation(self, intervention):
        info = intervention.demo_info
        self.intervention_demos[intervention] = []
        for d in self.demographics:
            size = info[d.name]["capacity"]
            vdemo = self.Demographic(
                d.name + "_V",
                size,
                size,
                d.cfr * info[d.name]["cfr_delta"],
                self.dur,
                intervention.quarantine,
                d,
                self.day,
            )
            d.children.append(vdemo)
            self.intervention_demos[intervention].append(vdemo)
            if intervention.quarantine:
                d.size -= size
                self.pop -= size

            for res in self.resources:
                res.demo_info[d.name + "_V"] = res.demo_info[d.name]

    def init_quarantine(self, intervention):
        info = intervention.demo_info
        self.intervention_demos[intervention] = []
        for d in self.demographics:
            qdemo = self.Demographic(
                d.name + "_Q",
                0,
                0,
                d.cfr * info[d.name]["cfr_delta"],
                self.dur,
                True,
                d,
                self.day,
                d.r_inf * info[d.name]["r_inf_delta"],
                d.r_sick * info[d.name]["r_sick_delta"],
            )
            d.children.append(qdemo)
            self.intervention_demos[intervention].append(qdemo)

    def end_intervention(self, intervention):
        if intervention.quarantine:
            for d in self.intervention_demos[intervention]:
                demo = d.parent
                demo.size += d.size
                self.pop += d.size
                demo.infected += d.infected
                self.infected += d.infected
                demo.sick += d.sick
                demo.immune += d.immune
                for i in range(len(d.case_hist)):
                    demo.case_hist[i] = (
                        demo.case_hist[i][0] + d.case_hist[i][0],
                        demo.case_hist[i][1],
                    )
                for i in range(len(d.inf_hist)):
                    demo.inf_hist[i] += d.inf_hist[i]
                    self.inf_hist[i] += d.inf_hist[i]
                    demo.death_hist[i] += d.death_hist[i]

            self.intervention_demos[intervention] = []

    def advance(self):
        self.day += 1
        for i in self.interventions:
            if i.crit == "Infected":
                if self.infected > i.threshold:
                    if not i.active:
                        if i.variolation:
                            self.init_variolation(i)
                        elif i.quarantine:
                            self.init_quarantine(i)
                    i.active = True
                else:
                    if i.active:
                        self.end_intervention(i)
                    i.active = False

            elif i.crit == "Day":
                if self.day in i.threshold:
                    if i.active:
                        i.active = False
                        self.end_intervention(i)
                    else:
                        i.active = True
                        if i.variolation:
                            self.init_variolation(i)
                        elif i.quarantine:
                            self.init_quarantine(i)

        dems = self.demographics
        if self.intervention_demos.keys():
            dems = list(
                set(dems).union(
                    demo
                    for demo in self.intervention_demos[i]
                    for i in self.intervention_demos.keys()
                )
            )
        for d in dems:
            self.adv_demo(d)

        self.inf_hist.append(self.infected)
        self.death_hist.append(self.dead)

    def adv_demo(self, d):
        r = self.r0
        cfr = d.cfr
        for res in self.resources:
            if res.used > res.n:
                cfr *= res.demo_info[d.name]["cfr_delta"]
        for i in self.interventions:
            if i.active and not i.variolation and not i.quarantine:
                r *= i.demo_info[d.name]["r_delta"]
                cfr *= i.demo_info[d.name]["cfr_delta"]

        concluding = d.case_hist[self.day]
        deaths = int(concluding[0] * concluding[1] * (1 - self.asymp))
        d.dead += deaths
        self.infected -= deaths
        d.infected -= deaths
        d.sick -= deaths
        self.dead += deaths
        d.death_hist.append(d.dead)

        new_sick = d.case_hist[-self.inc][0] * (1 - self.asymp)
        d.sick += new_sick

        new_cases = int(
            (d.size - d.immune - d.infected - d.dead) * (r * self.infected / self.dur) / self.pop
        )
        self.infected += new_cases
        d.infected += new_cases
        d.case_hist.append([new_cases, cfr])
        d.inf_hist.append(d.infected)

        recoveries = concluding[0] - deaths
        d.immune += recoveries
        self.infected -= recoveries
        d.infected -= recoveries
        d.sick -= recoveries * (1 - self.asymp)

        for res in self.resources:
            res.used += res.demo_info[d.name]["utilization"] * (
                new_sick - concluding[0] * (1 - self.asymp)
            )

    def adv_child_demo(self, d):
        r = self.r0
        cfr = d.cfr
        for res in self.resources:
            if res.used > res.n:
                cfr *= res.demo_info[d.name]["cfr_delta"]
        for i in self.interventions:
            if i.active and not i.variolation and not i.quarantine:
                r *= i.demo_info[d.name]["r_delta"]
                cfr *= i.demo_info[d.name]["cfr_delta"]

        concluding = d.case_hist[self.day]
        deaths = int(concluding[0] * concluding[1] * (1 - self.asymp))
        d.parent.dead += deaths
        d.parent.size += deaths
        d.size -= deaths
        if d.quarantine:
            self.pop += deaths
        else:
            self.infected -= deaths
        d.infected -= deaths
        d.sick -= deaths
        self.dead += deaths
        d.death_hist.append(d.dead)

        new_sick = d.case_hist[-self.inc][0] * (1 - self.asymp)
        d.sick += new_sick

        recoveries = concluding[0] - deaths
        d.parent.immune += recoveries
        d.parent.size += recoveries
        d.size -= recoveries
        if not d.quarantine:
            self.infected -= recoveries
        d.infected -= recoveries
        d.sick -= recoveries * (1 - self.asymp)

        new_cases = int(
            (d.size - d.immune - d.infected - d.dead) * (r * self.infected / self.dur) / self.pop
        )
        if d.name[-2:] == "_V":
            vulnerable = d.parent.size - d.parent.immune - d.parent.infected - d.parent.dead
            new_cases = min(concluding[0], vulnerable)
            d.size += new_cases
            d.parent.size -= new_cases
            if vulnerable < concluding[0]:
                self.pop += concluding[0] - vulnerable
        elif d.quarantine:
            space = d.cap - d.size
            if space > 0:
                admitting = min(d.parent.sick, space)
                d.parent.size -= admitting
                d.parent.infected -= admitting
                d.parent.sick -= admitting
                d.size += admitting
                d.infected += admitting
                d.sick += admitting
        if not d.quarantine:
            self.infected += new_cases
        d.infected += new_cases
        d.case_hist.append([new_cases, cfr])
        d.inf_hist.append(d.infected)

        for res in self.resources:
            res.used += res.demo_info[d.name]["utilization"] * (
                new_sick - concluding[0] * (1 - self.asymp)
            )

    def get_hist(self):
        data = {d.name: (d.inf_hist, d.death_hist) for d in self.demographics}
        data["Total"] = (self.inf_hist, self.death_hist)
        return data

    class Demographic:
        def __init__(
            self,
            name,
            size,
            infected,
            cfr,
            duration,
            quarantine=False,
            parent=None,
            day=0,
            r_inf=0,
            r_sick=0,
        ):
            self.name = name
            self.size = size
            self.infected = infected
            self.case_hist = [(0, 0) for i in range(duration - 1)]
            for _ in range(day - 1):
                self.case_hist.append((0, 0))
            self.case_hist.append((infected, cfr))
            self.sick = 0
            self.dead = 0
            self.immune = 0
            self.cfr = cfr
            self.r_inf = r_inf
            self.r_sick = r_sick
            self.quarantine = quarantine
            self.parent = parent
            self.children = []
            self.inf_hist = []
            self.death_hist = []

    class Intervention:
        # demo_info must contain r_delta, cfr_delta, and quarantine capacity for each demo
        def __init__(
            self, name, criterion, threshold, demo_info, variolation=False, quarantine=False
        ):
            self.name = name
            self.crit = criterion
            self.threshold = int(threshold)
            self.demo_info = {}
            for dem in demo_info.keys():
                self.demo_info[dem] = {}
                for key in demo_info[dem].keys():
                    self.demo_info[dem][key] = float(demo_info[dem][key])
            self.quarantine = quarantine
            self.variolation = variolation
            self.active = False

    class Resource:
        # demo_info is a dictionary of demographics, each containing a dictionary with utilization percentage and cfr_delta
        def __init__(self, name, number, demo_info):
            self.name = name
            self.n = int(number)
            self.used = 0
            self.demo_info = {}
            for dem in demo_info.keys():
                self.demo_info[dem] = {}
                for key in demo_info[dem].keys():
                    self.demo_info[dem][key] = float(demo_info[dem][key])


# %%

if __name__ == "__main__":
    inters = [
        {
            "name": "variolation",
            "criterion": "Infected",
            "threshold": 100,
            "demo_info": {
                "UnderSixty": {"capacity": 1000000, "cfr_delta": 0.3},
                "Sixties": {"capacity": 0, "cfr_delta": 1},
                "Seventies": {"capacity": 0, "cfr_delta": 1},
                "OverEighty": {"capacity": 0, "cfr_delta": 1},
            },
            "variolation": True,
            "quarantine": True,
        }
    ]
    demos = [
        {"name": "UnderSixty", "size": 235624706, "infected": 10000, "cfr": 0.01},
        {"name": "Sixties", "size": 20000940, "infected": 10000, "cfr": 0.036},
        {"name": "Seventies", "size": 15376083, "infected": 10000, "cfr": 0.08},
        {"name": "OverEighty", "size": 10420177, "infected": 10000, "cfr": 0.148},
    ]
    res = [
        {
            "name": "ventilators",
            "number": 62188,
            "demo_info": {
                "UnderSixty": {"utilization": 0.01, "cfr_delta": 3},
                "Sixties": {"utilization": 0.036, "cfr_delta": 3},
                "Seventies": {"utilization": 0.08, "cfr_delta": 3},
                "OverEighty": {"utilization": 0.148, "cfr_delta": 3},
            },
        }
    ]
    pop = Population(2.1, 10, 25, 0.5, demos, inters, res)
    for i in range(300):
        pop.advance()
    for i in pop.interventions:
        if i.active:
            pop.end_intervention(i)


# %%
