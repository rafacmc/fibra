import pandas as pd
import numpy as np

from datetime import date, timedelta, datetime
from pandas.tseries.offsets import DateOffset
from scipy.interpolate import CubicSpline
from bizdays import Calendar

cal = Calendar.load("../docs/ANBIMA.cal")

def switch(case):
            
    switcher = {
        0: "Zero",
        1: "Annual",
        2: "Semi-annual",
        3: "Fourth month",
        4: "Quarterly",
        6: "Bimonthly",
        12:"Monthly"
    }
        
    return switcher.get(case, "Invalid Argument")


class Bond:
    
    def __init__(self, date, maturity, ytm, coupon, freq):
        
        self.str_date = date
        self.str_maturity = maturity
        self.date = datetime.strptime(date, '%Y-%m-%d')
        self.maturity = datetime.strptime(maturity, '%Y-%m-%d')
        self.ytm = ytm
        self.coupon = coupon
        self.freq = freq
        self.frequency = switch(freq)
        self.convention = "BUS/252"
        self.nominal_value = 100
        self.coupon_rate = round(((1+coupon/100)**(1/self.freq)-1), 6)
        self.coupon_value = round(self.coupon_rate*self.nominal_value, 6)
        
        # Schedule
        period = int(12 / self.freq)
        initial_coupon = datetime(self.date.year-1, self.maturity.month, self.maturity.day)
        months = (self.maturity.year - initial_coupon.year)*12 + self.maturity.month - initial_coupon.month
        self.coupon_dates = pd.DatetimeIndex([self.maturity - DateOffset(months=e) for e in range(0, months, period)][::-1]).insert(0, initial_coupon).to_pydatetime()
        
        self.bond_schedule_official = list(filter(lambda x: x >= self.date, self.coupon_dates))
        self.bond_schedule = list([cal.adjust_next(d) for d in self.bond_schedule_official])
        
        # Cash flow
        self.du = [cal.bizdays(self.date, dt) for dt in self.bond_schedule]
        self.years = [float(f"{d/252:.14f}") for d in self.du]
        
        self.cfs = [self.coupon_value for d in self.du[:-1]] + [self.coupon_value + self.nominal_value]
        
        self.pvcfs = [round((self.coupon_value)/((1+self.ytm/100)**(float(f"{d/252:.10f}"))), 10) for d in self.du[:-1]] + [round((self.coupon_value+self.nominal_value)/((1+self.ytm/100)**(float(f"{self.du[-1]/252:.14f}"))), 10)]
        
        # Coupon dates
        self.last_coupon_date = list(filter(lambda x: x < self.date, self.coupon_dates))[-1]
        self.next_coupon_date = self.bond_schedule[0]
        self.days_to_next_coupon = self.du[0]
        self.num_payments = len(self.bond_schedule)
        
        # Term to maturity
        self.term_to_maturity = round(self.du[-1] / 252, 2)
        
        # Bond Price/Quote
        self.price = sum(self.pvcfs)
        
        # Current yield
        self.current_yield = round(((self.coupon / 100 * self.nominal_value) / self.price) * 100, 4)
        
        # Duration
        self.duration = sum([du * cf / self.price for du, cf in zip(self.du, self.pvcfs)])
        self.annual_duration = self.duration / 252
        self.modified_duration = round(self.annual_duration / (1 + self.ytm / 100), 4)
        
        # DV01 - Dollar duration
        self.dv01 = self.price * self.modified_duration / 10000
        
        # Convexity
        self.convexity = sum([(yrs + 1) * (yrs * cf / self.price) for yrs, cf in zip(self.years, self.pvcfs)]) * 1 / (1 + self.ytm / 100)**2
        
        # Price sensibility
        dy = 0.01
        self.sensibility = -self.modified_duration * dy + 1/2 * self.convexity * dy**2       
        
    
    def schedule(self):
                    
        return pd.DataFrame({"Date": self.bond_schedule})
    
    
    def schedule_official(self):
                
        return pd.DataFrame({"Date": self.bond_schedule_official})
        
    
    def cashflow(self):

        return pd.DataFrame({"Date": self.bond_schedule, 
                             "Days": self.du,
                             "Years": self.years,
                             "Cash Flows": self.cfs,
                             "PV of Cash Flows": self.pvcfs})
    
    
    def des(self):
        
        self.dict_bond = {"Convention": self.convention,
                          "Frequency": self.frequency,
                          "Current date": self.date.date(),
                          "Maturity date": self.maturity.date(),
                          "Last coupon date": self.last_coupon_date.date(),
                          "Next coupon date": self.next_coupon_date,
                          "Days to next coupon": self.days_to_next_coupon,
                          "Remaining payments": self.num_payments,
                          "Term to maturity": self.term_to_maturity,
                          "Coupon rate": self.coupon,
                          "Coupon value": self.coupon_value,
                          "Nominal value": self.nominal_value,
                          "Price": self.price,
                          "Yield to matutity": self.ytm,
                          "Current yield": self.current_yield,
                          "Duration": round(self.annual_duration, 4),
                          "Modified Duration": round(self.modified_duration, 4),
                          "DV01": round(self.dv01, 4),
                          "Convexity": round(self.convexity, 4),
                          "Sensibility": round(self.sensibility*100, 4)}
        
        return pd.DataFrame({"Bond":self.dict_bond}, index=self.dict_bond.keys())
        
        
    def get(self, item="price"):
        
        item_value = str(item).lower()
            
        self.bond_data = {"convention": self.convention,
                          "frequency": self.frequency,
                          "current_date": self.date.date(),
                          "maturity_date": self.maturity.date(),
                          "last_coupon_date": self.last_coupon_date.date(),
                          "next_coupon_date": self.next_coupon_date,
                          "days_to_next_coupon": self.days_to_next_coupon,
                          "remaining_payments": self.num_payments,
                          "term_to_maturity_years": self.term_to_maturity,
                          "coupon_rate": self.coupon,
                          "coupon_value": self.coupon_value,
                          "nominal_value": self.nominal_value,
                          "price": self.price,
                          "ytm": self.ytm,
                          "current_yield": self.current_yield,
                          "duration": self.annual_duration,
                          "modified_duration": self.modified_duration,
                          "dv01": self.dv01,
                          "convexity": self.convexity,
                          "sensibility": self.sensibility}
       
        return self.bond_data.get(item_value, "Invalid Item. Tap 'help(bond)' for more info.")
    
    
    def bond_curve(self):
        
        self.curve_dates = [dt for dt in cal.seq(self.date, self.maturity)]
        
        self.price_curve = [Bond(dt.strftime("%Y-%m-%d"), self.str_maturity, self.ytm, self.coupon, self.freq).get("price")                             for dt in self.curve_dates]
        
        series= pd.DataFrame({"Date": self.curve_dates, "Price": self.price_curve})
        
        return series
        

class ZeroBond:
    
    def __init__(self, date, maturity, ytm, nominal_value = 100):
        
        self.date = datetime.strptime(date, '%Y-%m-%d')
        self.maturity = datetime.strptime(maturity, '%Y-%m-%d')
        self.ytm = ytm
        self.convention = "BUS/252"
        self.nominal_value = nominal_value
        self.du = cal.bizdays(self.date, self.maturity)
    
        self.price = self.nominal_value / (1 + self.ytm / 100)**(self.du / 252)
        
        # Duration
        self.years = self.du / 252
        self.duration = self.years / (1 + self.ytm / 100)**(self.years + 1)
        self.annual_duration = self.duration / 252
        self.modified_duration = round(self.duration / (1 + self.ytm / 100), 4)
        
        # DV01 - Dollar duration
        self.dv01 = self.price * self.modified_duration / 10000
        
        # Convexity
        self.convexity = (self.years**2 + self.years) / (1 + self.ytm / 100)**(self.years+2)
        
        # Price sensibility
        dy = 0.01 # 1 bps
        self.sensibility = -self.modified_duration * dy + 1/2 * self.convexity * dy**2
    
    
    def des(self):
        
        self.dict_zero_bond = {"Convention": self.convention,
                               "Current date": self.date.date(),
                               "Maturity date": self.maturity.date(),
                               "Term to maturity": self.du,
                               "Nominal value": self.nominal_value,
                               "Price": self.price,
                               "Yield to maturity": self.ytm,
                               "Duration": round(self.duration, 4),
                               "Modified duration": round(self.modified_duration, 4),
                               "DV01": round(self.dv01, 4),
                               "Convexity": round(self.convexity, 4),
                               "Sensibility": round(self.sensibility*100, 4)}
        
        return pd.DataFrame({"Bond":self.dict_zero_bond}, index=self.dict_zero_bond.keys())
        
    
    def get(self, item = "price"):
        
        item_value = str(item).lower()
            
        self.bond_data = {"convention": self.convention,
                          "current_date": self.date.date(),
                          "maturity_date": self.maturity.date(),
                          "term_to_maturity_days": self.du,
                          "nominal_value": self.nominal_value,
                          "price": self.price,
                          "ytm": self.ytm,
                          "duration": self.duration,
                          "modified_duration": self.modified_duration,
                          "dv01": self.dv01,
                          "convexity": self.convexity,
                          "sensibility": self.sensibility}
       
        return self.bond_data.get(item_value, "Invalid Item. Tap 'help(bond)' for more info.")
    
    
    def bond_curve(self):
        
        self.curve_dates = [dt for dt in cal.seq(self.date, self.maturity)]
        
        self.price_curve = [ZeroBond(datetime.strftime(dt, "%Y-%m-%d"), datetime.strftime(self.maturity, "%Y-%m-%d"),                             self.ytm, self.nominal_value).get("price") for dt in self.curve_dates]
        
        series= pd.DataFrame({"Date": self.curve_dates, "Price": self.price_curve})
        
        return series


class YieldCurve:
    
    def __init__(self, vertices=None, rates=None, pos_date=date.today(), file_path=None, sep=","):
        
        self.date = pos_date
        self.vertices = vertices
        self.rates = rates
        self.file_path = str(file_path)
        self.sep = str(sep)

        thous, decim = (".", ",") if self.sep == ";" else (",", ".")
        
        if self.vertices is not None and self.rates is not None:
            if isinstance(self.vertices, list) and isinstance(self.rates, list):
                self.di_curve = pd.DataFrame(zip(self.vertices, self.rates), columns=["NSQ", "ULT"])
            else:
                raise TypeError("Parameters must be a list.")
        else:
            try:
                self.di_curve = pd.read_csv(str(self.file_path), sep=self.sep, thousands=thous, decimal=decim)
                self.di_curve.columns = ["NSQ", "ULT"]
            except:
                raise TypeError("Specified file path does not exist or file type is not supported.")
        
    
    def load(self):
        
        return self.di_curve
    
    
    def curve(self, method="cs"):
        
        self.di_vert = self.di_curve["NSQ"].values
        self.di_rate = self.di_curve["ULT"].values
        self.interpolation_method = method
        
        self.daily_sequence = [datetime.strftime(dt, "%Y-%m-%d") for dt in cal.seq(self.date, date(2035, 1, 1))]
        
        xs = [cal.bizdays(self.date, dt) for dt in self.daily_sequence]
        
        if self.interpolation_method == "cs" or self.interpolation_method == "cubic_spline":
            cs = CubicSpline(self.di_vert, self.di_rate)
            self.interpolate_rate = cs(xs)
        else:
            TypeError("Invalid interpolation method.")
        
        self.di_curve_interpolate = pd.DataFrame({"Date": self.daily_sequence, "DI Futures": self.interpolate_rate})

        self.du_curve = [cal.bizdays(self.date, dt) for dt in self.di_curve_interpolate["Date"].values]
        self.di_curve_interpolate.insert(1, "Days", self.du_curve)

        self.di_factor = [((1 + rate / 100)**(du/252)) for rate, du in zip(self.interpolate_rate, self.du_curve)]
        self.di_curve_interpolate.insert(2, "Factor", self.di_factor)

        self.di_curve_interpolate.set_index("Date", inplace=True)
        
        return self.di_curve_interpolate
    
    
    def forward(self, start_date, end_date):
        
        self.start_date = cal.adjust_next(start_date)
        self.end_date = cal.adjust_next(end_date)
        
        self.start_date = datetime.strftime(self.start_date, '%Y-%m-%d')
        self.end_date = datetime.strftime(self.end_date, '%Y-%m-%d')
        
        self.start_point_rate = self.di_curve_interpolate.loc[self.start_date, "Factor"]
        self.end_point_rate = self.di_curve_interpolate.loc[self.end_date, "Factor"]
        
        self.start_point_du = self.di_curve_interpolate.loc[self.start_date, "Days"]
        self.end_point_du = self.di_curve_interpolate.loc[self.end_date, "Days"]
        self.term = self.end_point_du - self.start_point_du
        
        self.forward_rate = ((self.end_point_rate / self.start_point_rate)**(252 / self.term) - 1) * 100
        
        return round(self.forward_rate, 8)


class DebentureDI:
    
    def __init__(self, date, maturity, vne, vna, pu, issue_rate, market_rate, freq, redemption, payment="adjust", vertices=None, exp_rates=None, yield_curve_file=None, sep=","):
        
        self.date = datetime.strptime(date, '%Y-%m-%d')
        self.date_str = date.replace("-", "")
        self.maturity = datetime.strptime(maturity, '%Y-%m-%d')
        self.vne = vne
        self.vna = vna
        self.pu = pu
        self.issue_rate = issue_rate
        self.market_rate = market_rate
        self.freq = freq
        self.redempt = redemption
        self.payment = payment
        self.vertices = vertices
        self.rates = exp_rates
        self.file_path = yield_curve_file
        
        # Schedule
        period = int(12 / self.freq)
        initial_event = datetime(self.date.year-1, self.maturity.month, self.maturity.day)
        months = (self.maturity.year - initial_event.year)*12 + self.maturity.month - initial_event.month
        self.event_dates = pd.DatetimeIndex([self.maturity - DateOffset(months=e) for e in range(0, months, period)][::-1]).insert(0, initial_event).to_pydatetime()
        
        self.debent_schedule_official = list(filter(lambda x: x >= self.date, self.event_dates))
        self.debent_schedule = list([cal.adjust_next(d) for d in self.debent_schedule_official])
        self.str_date_schedule = [datetime.strftime(dt, "%Y-%m-%d") for dt in self.debent_schedule]
        
        # Create a worday calendar
        self.du = [cal.bizdays(self.date, dt) for dt in self.debent_schedule]
        self.inter_du = [cal.bizdays(dt1, dt) for dt1, dt in zip(self.debent_schedule, self.debent_schedule[1:])]
        
        # Redemption rate
        self.df_sechedule = pd.DataFrame(self.str_date_schedule, columns=["Date"])
        self.df_redempt_schedule = pd.DataFrame(list(self.redempt.items()), columns=["Date", "Redemption"])
        self.df_redempt = pd.merge(self.df_sechedule, self.df_redempt_schedule, on="Date", how="left").fillna(0)
        
        self.redempt_per = self.df_redempt["Redemption"].values
        
        # Adjust the VNA
        self.vna_val = [self.vna] * len(self.redempt_per)
        self.vna_adj = [self.vna]
        
        if self.payment == "fixed":
            [self.vna_adj.append(self.vna_adj[i] - self.vna_adj[0]*r) for i, r in zip(range(0, len(self.redempt_per)), self.redempt_per)]
            self.vna_adj = self.vna_adj[:-1]
        else:
            [self.vna_adj.append(self.vna_adj[i]*(1-r)) for i, r in zip(range(0, len(self.redempt_per)), self.redempt_per)]
            self.vna_adj = self.vna_adj[:-1]
        
        # Redemption value
        if self.payment == "fixed":
            self.redemption_val = [self.vna*redempt for redempt in self.redempt_per]
        else:
            self.redemption_val = [vna_adj*redempt for vna_adj, redempt in zip(self.vna_adj, self.redempt_per)]
        
        # DI projections
        self.di_interpolate = YieldCurve(pos_date=date, vertices=self.vertices, rates=self.rates, file_path=self.file_path, sep=";").curve()
        self.di_proj = [self.di_interpolate.loc[datetime.strftime(dt, "%Y-%m-%d"), "DI Futures"] for dt in self.debent_schedule]
        
        # Calculate the factor
        self.factor = [(((1 + i / 100)**(1 / 252) - 1)*self.issue_rate / 100 + 1 )**d for i, d in zip(self.di_proj, self.du)]
        
        # Calculate the term
        self.f_term = ((((1 + self.di_proj[0] / 100)**(1 / 252) - 1)*self.issue_rate / 100 + 1)*((((1 + self.di_proj[0] / 100)**(1 / 252) - 1)*self.issue_rate / 100 + 1)**self.du[0]) - 1)
        
        self.term = [self.f_term] + [(f1 / f0) - 1 for f1, f0 in zip(self.factor[1:], self.factor)]
        
        # Cash Flow
        self.interest = [vna * t for vna, t in zip(self.vna_adj, self.term)]
        self.cfs = [j + a for j, a in zip(self.interest, self.redemption_val)]
        
        # Discount factor
        self.discount = [1 / (((1 + i / 100)**(1 / 252) - 1)*self.market_rate / 100 + 1 )**d for i, d in zip(self.di_proj, self.du)]
        
        # PV of Cash flows
        self.pvcfs = [cf*discount for cf, discount in zip(self.cfs, self.discount)]
        
        # Price
        self.price = sum(self.pvcfs)
        
        # Price ratio
        self.price_pu_ratio = self.price / self.pu
        
        # Payments flow
        self.last_payment_date = list(filter(lambda x: x < self.date, self.event_dates))[-1]
        self.next_payment_date = self.debent_schedule[0]
        self.days_to_next_payment = self.du[0]
        self.num_payments = len(self.debent_schedule)

        # Term to maturity
        self.term_to_maturity = round(self.du[-1] / 252, 2)
        
        # Duration
        self.duration = sum([du * cf / self.price for du, cf in zip(self.du, self.pvcfs)])
        self.annual_duration = self.duration / 252
    
    
    def des(self):
        
        self.dict_debent = {"Convention": "Percent DI",
                            "Frequency": switch(self.freq),
                            "Date": self.date.date(),
                            "Maturity": self.maturity.date(),
                            "Last payment date": self.last_payment_date.date(),
                            "Next payment date": self.next_payment_date,
                            "Days to next payment": self.days_to_next_payment,
                            "Remaining payments": self.num_payments,
                            "Term to maturity": self.term_to_maturity,
                            "VNE": "{:.6f}".format(self.vne),
                            "VNA": "{:.6f}".format(self.vna),
                            "PU Par": "{:.6f}".format(self.pu),
                            "Price": "{:.6f}".format(self.price),
                            "Price/Par ratio": "{:.2%}".format(self.price_pu_ratio),
                            "Issue rate": "{:.4%}".format(self.issue_rate / 100),
                            "Market rate": "{:.4%}".format(self.market_rate / 100),
                            "Duration": round(self.annual_duration, 2)}
        
        return pd.DataFrame({"Debenture DI":self.dict_debent}, index=self.dict_debent.keys())        
        
    
    def cashflow(self):
        
        return pd.DataFrame({"Date": self.debent_schedule,
                             "Days": self.du,
                             "DI": self.di_proj,
                             "Interests": self.interest,
                             "Redemptions": self.redemption_val,
                             "Payments": self.cfs,
                             "PV of Cash Flow": self.pvcfs})

    
    def get(self, item="price"):
        
        item_value = str(item).lower()
            
        self.debent_data = {"convention": "Percent DI",
                            "frequency": self.freq,
                            "current_date": self.date.date(),
                            "maturity_date": self.maturity.date(),
                            "last_payment_date": self.last_payment_date.date(),
                            "next_payment_date": self.next_payment_date,
                            "days_to_next_payment": self.days_to_next_payment,
                            "remaining_payments": self.num_payments,
                            "term_to_maturity_years": self.term_to_maturity,
                            "vne": "{:.6f}".format(self.vne),
                            "vna": "{:.6f}".format(self.vna),
                            "pu_par": "{:.6f}".format(self.pu),
                            "price": "{:.6f}".format(self.price),
                            "issue_rate": "{:.4%}".format(self.issue_rate / 100),
                            "marker_rate": "{:.4%}".format(self.market_rate / 100),
                            "duration": round(self.annual_duration, 2)}
       
        return self.debent_data.get(item_value, "Invalid Item. Tap 'help(debent)' for more info.")


class DebentureSpread:
    
    def __init__(self, date, maturity, vne, vna, pu, issue_spread, market_spread, freq, redemption, payment="adjust", vertices=None, exp_rates=None, yield_curve_file=None, sep=","):
        
        self.date = datetime.strptime(date, '%Y-%m-%d')
        self.date_str = date.replace("-", "")
        self.maturity = datetime.strptime(maturity, '%Y-%m-%d')
        self.vne = vne
        self.vna = vna
        self.pu = pu
        self.issue_spread = issue_spread
        self.market_spread = market_spread
        self.freq = freq
        self.redempt = redemption
        self.payment = payment
        self.vertices = vertices
        self.rates = exp_rates
        self.file_path = yield_curve_file
        
        # Schedule
        period = int(12 / self.freq)
        initial_event = datetime(self.date.year-1, self.maturity.month, self.maturity.day)
        months = (self.maturity.year - initial_event.year)*12 + self.maturity.month - initial_event.month
        self.event_dates = pd.DatetimeIndex([self.maturity - DateOffset(months=e) for e in range(0, months, period)][::-1]).insert(0, initial_event).to_pydatetime()
        
        self.debent_schedule_official = list(filter(lambda x: x >= self.date, self.event_dates))
        self.debent_schedule = list([cal.adjust_next(d) for d in self.debent_schedule_official])
        self.str_date_schedule = [datetime.strftime(dt, "%Y-%m-%d") for dt in self.debent_schedule]
        self.first_du = cal.bizdays(self.date, self.debent_schedule[0])
        
        # Create a workday calendar
        self.du = [cal.bizdays(self.date, dt) for dt in self.debent_schedule]
        self.inter_du = [cal.bizdays(dt1, dt) for dt1, dt in zip(self.debent_schedule, self.debent_schedule[1:])]
        
        # Redemption rate
        self.df_sechedule = pd.DataFrame(self.str_date_schedule, columns=["Date"])
        self.df_redempt_schedule = pd.DataFrame(list(self.redempt.items()), columns=["Date", "Redemption"])
        self.df_redempt = pd.merge(self.df_sechedule, self.df_redempt_schedule, on="Date", how="left").fillna(0)
        
        self.redempt_per = self.df_redempt["Redemption"].values
        
        # Adjust the VNA
        self.vna_val = [self.vna] * len(self.redempt_per)
        self.vna_adj = [self.vna]
        
        if self.payment == "fixed":
            [self.vna_adj.append(self.vna_adj[i] - self.vna_adj[0]*r) for i, r in zip(range(0, len(self.redempt_per)), self.redempt_per)]
            self.vna_adj = self.vna_adj[:-1]
        else:
            [self.vna_adj.append(self.vna_adj[i]*(1-r)) for i, r in zip(range(0, len(self.redempt_per)), self.redempt_per)]
            self.vna_adj = self.vna_adj[:-1]
        
        # Redemption value
        if self.payment == "fixed":
            self.redemption_val = [self.vna*redempt for redempt in self.redempt_per]
        else:
            self.redemption_val = [vna_adj*redempt for vna_adj, redempt in zip(self.vna_adj, self.redempt_per)]
        
        # DI projections
        self.di_interpolate = YieldCurve(pos_date=date, vertices=self.vertices, rates=self.rates, file_path=self.file_path, sep=";").curve()
        self.di_proj = [self.di_interpolate.loc[datetime.strftime(dt, "%Y-%m-%d"), "DI Futures"] for dt in self.debent_schedule]
            
        # Calculate the factor
        self.f_factor = self.pu*(((((1 + self.di_proj[0] / 100)**(1 / 252))*((1 + self.issue_spread / 100)**(1 / 252)))**self.first_du) - 1)
        self.factor = [(((1 + i / 100)**(1 / 252))*((1 + self.issue_spread / 100)**(1 / 252)))**d for i, d in zip(self.di_proj, self.du)]
            
        self.term = [(f1 / f0) - 1 for f1, f0 in zip(self.factor[1:], self.factor)]
        
        # Cash Flow
        self.interest = [(self.pu - self.vna) + self.f_factor] + [vna * t for vna, t in zip(self.vna_adj, self.term)]
        self.cfs = [j + a for j, a in zip(self.interest, self.redemption_val)]
        
        # Discount factor
        self.discount = [1 / (((1 + i / 100)**(1 / 252))*((1 + self.market_spread / 100)**(1 / 252)))**d for i, d in zip(self.di_proj, self.du)]
        
        # PV of Cash flows
        self.pvcfs = [cf*discount for cf, discount in zip(self.cfs, self.discount)]
        
        # Price
        self.price = sum(self.pvcfs)
        
        # Price ratio
        self.price_pu_ratio = self.price / self.pu
        
        # Payments flow
        self.last_payment_date = list(filter(lambda x: x < self.date, self.event_dates))[-1]
        self.next_payment_date = self.debent_schedule[0]
        self.days_to_next_payment = self.du[0]
        self.num_payments = len(self.debent_schedule)

        # Term to maturity
        self.term_to_maturity = round(self.du[-1] / 252, 2)
        
        # Duration
        self.duration = sum([du * cf / self.price for du, cf in zip(self.du, self.pvcfs)])
        self.annual_duration = self.duration / 252
    
    
    def des(self):
        
        self.dict_debent = {"Convention": "DI Spread",
                            "Frequency": switch(self.freq),
                            "Date": self.date.date(),
                            "Maturity": self.maturity.date(),
                            "Last payment date": self.last_payment_date.date(),
                            "Next payment date": self.next_payment_date,
                            "Days to next payment": self.days_to_next_payment,
                            "Remaining payments": self.num_payments,
                            "Term to maturity": self.term_to_maturity,
                            "VNE": "{:.6f}".format(self.vne),
                            "VNA": "{:.6f}".format(self.vna),
                            "PU Par": "{:.6f}".format(self.pu),
                            "Price": "{:.6f}".format(self.price),
                            "Price/Par ratio": "{:.2%}".format(self.price_pu_ratio),
                            "Issue rate": "{:.4%}".format(self.issue_spread / 100),
                            "Market rate": "{:.4%}".format(self.market_spread / 100),
                            "Duration": round(self.annual_duration, 2)}
        
        return pd.DataFrame({"Debenture DI Spread":self.dict_debent}, index=self.dict_debent.keys())        
        
    
    def cashflow(self):
        
        return pd.DataFrame({"Date": self.debent_schedule,
                             "Days": self.du,
                             "DI": self.di_proj,
                             "VNA": self.vna_adj,
                             "Interests": self.interest,
                             "Redemptions": self.redemption_val,
                             "Payments": self.cfs,
                             "PV of Cash Flows": self.pvcfs})

    
    def get(self, item="price"):
        
        item_value = str(item).lower()
            
        self.debent_data = {"convention": "DI Spread",
                            "frequency": self.freq,
                            "current_date": self.date.date(),
                            "maturity_date": self.maturity.date(),
                            "last_payment_date": self.last_payment_date.date(),
                            "next_payment_date": self.next_payment_date,
                            "days_to_next_payment": self.days_to_next_payment,
                            "remaining_payments": self.num_payments,
                            "term_to_maturity_years": self.term_to_maturity,
                            "vne": "{:.6f}".format(self.vne),
                            "vna": "{:.6f}".format(self.vna),
                            "pu_par": "{:.6f}".format(self.pu),
                            "price": "{:.6f}".format(self.price),
                            "issue_rate": "{:.4%}".format(self.issue_spread / 100),
                            "marker_rate": "{:.4%}".format(self.market_spread / 100),
                            "duration": round(self.annual_duration, 2)}
       
        return self.debent_data.get(item_value, "Invalid Item. Tap 'help(debent)' for more info.")


class DebentureIPCA:
    
    def __init__(self, date, maturity, vne, vna, pu, issue_rate, market_rate, freq, redemption, payment="adjust"):
        
        self.date = datetime.strptime(date, '%Y-%m-%d')
        self.maturity = datetime.strptime(maturity, '%Y-%m-%d')
        self.vne = vne
        self.vna = vna
        self.pu = pu
        self.issue_rate = issue_rate
        self.market_rate = market_rate
        self.freq = freq
        self.redempt = redemption
        self.payment = payment
        
        # Schedule
        period = int(12 / self.freq)
        initial_event = datetime(self.date.year-1, self.maturity.month, self.maturity.day)
        months = (self.maturity.year - initial_event.year)*12 + self.maturity.month - initial_event.month
        self.event_dates = pd.DatetimeIndex([self.maturity - DateOffset(months=e) for e in range(0, months, period)][::-1]).insert(0, initial_event).to_pydatetime()
        
        self.debent_schedule_official = list(filter(lambda x: x >= self.date, self.event_dates))
        self.debent_schedule = list([cal.adjust_next(d) for d in self.debent_schedule_official])
        self.str_date_schedule = [datetime.strftime(dt, "%Y-%m-%d") for dt in self.debent_schedule]
        self.first_du = cal.bizdays(self.date, self.debent_schedule[0])
        
        # Create a workday calendar
        self.du = [cal.bizdays(self.date, dt) for dt in self.debent_schedule]
        self.inter_du = [cal.bizdays(dt1, dt) for dt1, dt in zip(self.debent_schedule, self.debent_schedule[1:])]
        
        # Redemption rate
        self.df_sechedule = pd.DataFrame(self.str_date_schedule, columns=["Date"])
        self.df_redempt_schedule = pd.DataFrame(list(self.redempt.items()), columns=["Date", "Redemption"])
        self.df_redempt = pd.merge(self.df_sechedule, self.df_redempt_schedule, on="Date", how="left").fillna(0)
        
        self.redempt_per = self.df_redempt["Redemption"].values
        
        # Adjust the VNA
        self.vna_val = [self.vna] * len(self.redempt_per)
        self.vna_adj = [self.vna]
        
        if self.payment == "fixed":
            self.vna_adj = [self.vna] + [self.vna_val[i+1] - i*(self.vna_val[i]*self.redempt_per[i]) for i in range(len(self.redempt_per)-1)]
        else:
            [self.vna_adj.append(self.vna_adj[i]*(1-r)) for i, r in zip(range(0, len(self.redempt_per)), self.redempt_per)]
        
        # Shift VNA to calculate interest
        self.vna_for_interest = self.vna_adj[1:] + [self.vna_adj[-1]]
        
        # Redemption value
        if self.payment == "fixed":
            self.redemption_val = [self.vna*redempt for redempt in self.redempt_per]
        else:
            self.redemption_val = [vna_adj*redempt for vna_adj, redempt in zip(self.vna_adj, self.redempt_per)]
        
        # Calculate the interest
        self.f_factor = ((1 + self.issue_rate / 100)**(self.first_du / 252) - 1)
        self.factor = [((1 + self.issue_rate / 100)**(d / 252) - 1) for d in self.inter_du]
    
        self.f_interest = (self.pu - self.vna) + self.pu*self.f_factor
        self.interest = [self.f_interest] + [abs(vna_adj*f) for f, vna_adj in zip(self.factor, self.vna_for_interest)]
        
         # Cash flows
        self.cfs = [j + a for j, a in zip(self.interest, self.redemption_val)]
        
        # Discount factor
        self.discount = [1 / ((1 + self.market_rate / 100)**(d / 252)) for d in self.du]
        
        # PV of Cash flows
        self.pvcfs = [cf*discount for cf, discount in zip(self.cfs, self.discount)]
        
        # Price
        self.price = sum(self.pvcfs)
        
        # Price ratio
        self.price_pu_ratio = self.price / self.pu

        # Payments flow
        self.last_payment_date = list(filter(lambda x: x < self.date, self.event_dates))[-1]
        self.next_payment_date = self.debent_schedule[0]
        self.days_to_next_payment = self.du[0]
        self.num_payments = len(self.debent_schedule)

        # Term to maturity
        self.term_to_maturity = round(self.du[-1] / 252, 2)
        
        # Duration
        self.duration = sum([du * cf / self.price for du, cf in zip(self.du, self.pvcfs)])
        self.annual_duration = self.duration / 252
    
    
    def des(self):
        
        self.dict_debent = {"Convention": "IPCA Spread",
                            "Frequency": switch(self.freq),
                            "Date": self.date.date(),
                            "Maturity": self.maturity.date(),
                            "Last payment date": self.last_payment_date.date(),
                            "Next payment date": self.next_payment_date,
                            "Days to next payment": self.days_to_next_payment,
                            "Remaining payments": self.num_payments,
                            "Term to maturity": self.term_to_maturity,
                            "VNE": "{:.6f}".format(self.vne),
                            "VNA": "{:.6f}".format(self.vna),
                            "PU Par": "{:.6f}".format(self.pu),
                            "Price": "{:.6f}".format(self.price),
                            "Price/Par ratio": "{:.2%}".format(self.price_pu_ratio),
                            "Issue rate": "{:.4%}".format(self.issue_rate / 100),
                            "Market rate": "{:.4%}".format(self.market_rate / 100),
                            "Duration": round(self.annual_duration, 2)}
        
        return pd.DataFrame({"Debenture IPCA Spread":self.dict_debent}, index=self.dict_debent.keys())        
        
    
    def cashflow(self):
        
        return pd.DataFrame({"Date": self.debent_schedule,
                             "Days": self.du,
                             "Interests": self.interest,
                             "Redemptions": self.redemption_val,
                             "Payments": self.cfs,
                             "PV of Cash Flows": self.pvcfs})

    
    def get(self, item="price"):
        
        item_value = str(item).lower()
            
        self.debent_data = {"convention": "IPCA Spread",
                            "frequency": self.freq,
                            "current_date": self.date.date(),
                            "maturity_date": self.maturity.date(),
                            "last_payment_date": self.last_payment_date.date(),
                            "next_payment_date": self.next_payment_date,
                            "days_to_next_payment": self.days_to_next_payment,
                            "remaining_payments": self.num_payments,
                            "term_to_maturity_years": self.term_to_maturity,
                            "vne": "{:.6f}".format(self.vne),
                            "vna": "{:.6f}".format(self.vna),
                            "pu_par": "{:.6f}".format(self.pu),
                            "price": "{:.6f}".format(self.price),
                            "issue_rate": "{:.4%}".format(self.issue_rate / 100),
                            "marker_rate": "{:.4%}".format(self.market_rate / 100),
                            "duration": round(self.annual_duration, 2)}
       
        return self.debent_data.get(item_value, "Invalid Item. Tap 'help(debent)' for more info.")