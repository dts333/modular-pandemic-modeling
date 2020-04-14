""" 
Scarce resources are assigned as cases come in and cases come in in the order demographics are listed.
If a vulnerable group is more likely to be given the resources, they should be listed first.
"""
#%%

#%%
class Population:
    def __init__(
        self,
        r_inf,
        r_sick,
        incubation,
        duration,
        asymptomatic,
        demographics,
        interventions,
        resources,
    ):
        self.day = 0
        self.pop = 0
        self.infected = 0
        self.sick = 0
        self.dead = 0
        self.r_inf = r_inf
        self.r_sick = r_sick
        self.inc = incubation
        self.dur = duration
        self.sick_dur = duration - incubation
        self.asymp = asymptomatic
        self.demographics = []
        self.intervention_demos = {}
        self.interventions = []
        self.resources = []
        self.inf_hist = []
        self.death_hist = []

        for r in resources:
            if r["name"]:
                self.resources.append(self.Resource(r["name"], r["number"], r["demo_info"]))

        for d in demographics:
            self.pop += int(d["size"])
            self.infected += int(d["infected"])
            self.demographics.append(
                self.Demographic(
                    name=d["name"],
                    size=int(d["size"]),
                    infected=int(d["infected"]),
                    cfr=float(d["cfr"]),
                    duration=duration,
                    resources=self.resources,
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

    def advance(self):
        self.day += 1
        for i in self.interventions:
            if i.crit == "Infected":
                if self.infected > i.threshold:
                    if not i.active:
                        i.activate(self.demographics, self.resources, self.day, self.dur)
                else:
                    if i.active:
                        i.deactivate()

            elif i.crit == "Day":
                if self.day in i.threshold:
                    if i.active:
                        i.deactivate()
                    else:
                        i.activate(self.demographics, self.resources, self.day, self.dur)

        for d in self.demographics:
            self.adv_demo(d)

        self.inf_hist.append(self.infected)
        self.death_hist.append(self.dead)

    def adv_demo(self, d):
        cfr = d.cfr
        for child in d.children:
            self.adv_child_demo(child)

        for i in self.interventions:
            if i.active and not i.variolation and not i.quarantine:
                r *= i.demo_info[d.name]["r_delta"]
                cfr *= i.demo_info[d.name]["cfr_delta"]

        concluding = d.case_hist[self.day]
        deaths = int(concluding[0] * concluding[1] * (1 - self.asymp))
        d.dead += deaths
        self.infected -= deaths
        d.infected -= deaths
        self.sick -= deaths
        d.sick -= deaths
        self.dead += deaths
        d.death_hist.append(d.dead)

        new_sick = d.case_hist[-self.inc][0] * (1 - self.asymp)
        d.sick += new_sick
        self.sick += new_sick
        res_dict = {}
        for res in self.resources:
            res.used -= concluding[2][res.name]
            res_avail = res.n - res.used
            try:
                cfr = cfr * (
                    1
                    + (res.demo_info[d.name]["cfr_delta"] - 1)
                    * max(new_sick - res_avail, 0)
                    / new_sick
                )
            except ZeroDivisionError:
                pass
            deploying = min(new_sick * res.demo_info[d.name]["utilization"], res_avail)
            res_dict[res.name] = deploying
            res.used += deploying
        d.case_hist[-self.inc].append(cfr)
        d.case_hist[-self.inc].append(res_dict)

        new_cases = int(
            (d.size - d.immune - d.infected - d.dead)
            * (
                (self.r_inf * (self.infected - self.sick) / self.inc)
                + (self.r_sick * self.sick / self.sick_dur)
            )
            / self.pop
        )
        self.infected += new_cases
        d.infected += new_cases
        d.case_hist.append([new_cases])
        d.inf_hist.append(d.infected)

        recoveries = concluding[0] - deaths
        d.immune += recoveries
        self.infected -= recoveries
        d.infected -= recoveries
        self.sick -= recoveries * (1 - self.asymp)
        d.sick -= recoveries * (1 - self.asymp)

    def adv_child_demo(self, d):
        cfr = d.cfr
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
            self.sick -= deaths
        d.infected -= deaths
        d.sick -= deaths
        self.dead += deaths
        d.death_hist.append(d.dead)

        new_sick = d.case_hist[-self.inc][0] * (1 - self.asymp)
        d.sick += new_sick
        if not d.quarantine:
            self.sick += new_sick
        res_dict = {}
        for res in self.resources:
            res.used -= concluding[2][res.name]
            res_avail = res.n - res.used
            try:
                cfr = cfr * (
                    1
                    + (res.demo_info[d.name]["cfr_delta"] - 1)
                    * max(new_sick - res_avail, 0)
                    / new_sick
                )
            except ZeroDivisionError:
                pass
            deploying = min(new_sick * res.demo_info[d.name]["utilization"], res_avail)
            res_dict[res.name] = deploying
            res.used += deploying
        d.case_hist[-self.inc].append(cfr)
        d.case_hist[-self.inc].append(res_dict)

        recoveries = concluding[0] - deaths
        d.parent.immune += recoveries
        d.parent.size += recoveries
        d.size -= recoveries
        d.infected -= recoveries
        d.sick -= recoveries * (1 - self.asymp)
        if d.quarantine:
            self.pop += recoveries
        else:
            self.infected -= recoveries
            self.sick -= recoveries * (1 - self.asymp)

        new_cases = int(
            (d.size - d.immune - d.infected - d.dead)
            * ((self.r_inf * self.infected / self.inc) + (self.r_sick * self.sick / self.sick_dur))
            / self.pop
        )
        if d.name[-2:] == "_V":
            vulnerable = d.parent.size - d.parent.immune - d.parent.infected - d.parent.dead
            new_cases = min(concluding[0], vulnerable)
            d.size += new_cases
            d.parent.size -= new_cases
            if d.quarantine:
                self.pop -= new_cases
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
        d.case_hist.append([new_cases])
        d.inf_hist.append(d.infected)

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
            resources,
            quarantine=False,
            parent=None,
            cap=None,
            day=0,
            r_inf=0,
            r_sick=0,
        ):
            self.name = name
            self.size = size
            self.infected = infected
            res_dict = {r.name: 0 for r in resources}
            self.case_hist = [[0, 0, res_dict] for i in range(duration - 1)]
            for _ in range(day - 1):
                self.case_hist.append([0, 0, res_dict])
            self.case_hist.append([infected])
            self.sick = 0
            self.dead = 0
            self.immune = 0
            self.cfr = cfr
            self.r_inf = r_inf
            self.r_sick = r_sick
            self.quarantine = quarantine
            self.parent = parent
            self.cap = cap
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
            self.demos = []
            for dem in demo_info.keys():
                self.demo_info[dem] = {}
                for key in demo_info[dem].keys():
                    self.demo_info[dem][key] = float(demo_info[dem][key])
            self.quarantine = quarantine
            self.variolation = variolation
            self.active = False

        def activate(self, demos, resources, day, duration):
            self.active = True
            if self.quarantine:
                token = "_Q"
                size = 0
            elif self.variolation:
                token = "_V"
            for d in demos:
                if self.variolation:
                    size = self.demo_info[d.name]["capacity"]
                demo = self.Demographic(
                    name=d.name + token,
                    size=size,
                    infected=size,
                    cfr=d.cfr * self.demo_info[d.name]["cfr_delta"],
                    duration=duration,
                    resources=resources,
                    quarantine=self.quarantine,
                    parent=d,
                    cap=self.demo_info[d.name]["cap"],
                    day=day,
                    r_inf=d.r_inf * self.demo_info[d.name]["r_inf_delta"],
                    r_sick=d.r_sick * self.demo_info[d.name]["r_sick_delta"],
                )
                d.children.append(demo)
                self.demos.append(demo)

                if self.quarantine:
                    d.size -= size
                    self.pop -= size

                for res in resources:
                    res.demo_info[d.name + token] = res.demo_info[d.name]

        def deactivate(self):
            for d in self.demos:
                demo = d.parent
                demo.size += d.size
                demo.infected += d.infected
                demo.sick += d.sick
                demo.immune += d.immune
                if self.quarantine:
                    self.pop += d.size
                    self.infected += d.infected
                    self.sick += d.sick
                for i in range(len(d.case_hist)):
                    demo.case_hist[i][0] += d.case_hist[i][0]
                for i in range(len(d.inf_hist)):
                    demo.inf_hist[i] += d.inf_hist[i]
                    self.inf_hist[i] += d.inf_hist[i]
                    demo.death_hist[i] += d.death_hist[i]

            self.demos = []

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
