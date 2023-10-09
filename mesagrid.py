import os
import numpy as np
import pandas as pd
import scipy as sp
from tomso import gyre
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

class Track:
    def __init__(self, dir, parameters=None, 
                 load_history_extras=None, 
                 usecols_profiles=None,
                 usecols_history=None, 
                 cpus=None,
                 dir_name="LOGS"):
        if dir_name is not None:
            self.dir = os.path.join(dir, dir_name)
        else:
            self.dir = dir
        self.parameters = parameters
        self.cpus = os.cpu_count() if cpus is None else cpus
        
        self.usecols_profiles = usecols_profiles
        self.usecols_history = usecols_history
        
        self.load_history_extras = load_history_extras
        self.loaded = False
        
        self._history  = None
        self._index    = None
        self._profiles = None
        self._gyres    = None
        self._freqs    = None
    
    #def __repr__(self):
    #    return f"Track with {len(self.history)} models\nParameters: " + str(self.parameters)
    
    ### Lazy-load getters and setters for the data 
    @property
    def history(self):
        if self._history is None:
            self._history = self.load_history_file()
        return self._history
    
    @property
    def index(self):
        if self._index is None:
            self._index = self.get_index()
        return self._index
    
    # this should optionally take a set of profile numbers, and only load those profile numbers 
    # but it should still be accessible like a list??? like track.profiles[5] 
    @property
    def profiles(self):
        if self._profiles is None:
            self._profiles = self.get_profiles()
        return self._profiles
    
    # this should optionally take a set of profile numbers, and only load those profile numbers 
    @property
    def gyres(self):
        if self._gyres is None:
            self._gyres = self.get_gyres()
        return self._gyres
    
    @property
    def freqs(self):
        if self._freqs is None:
            self._freqs = self.get_freqs()
        return self._freqs
    
    ### HISTORY
    def load_history_file(self):
        print("Parsing", self.parameters)
        DF_ = pd.read_table(os.path.join(self.dir, 'history.data'), 
                            skiprows=5, sep='\s+', 
                            usecols=self.usecols_history)
        
        if self.load_history_extras is not None and not self.loaded:
            DF_ = self.load_history_extras(self, DF_)
        
        self.loaded = True
        return DF_
    
    ### INDEX
    def get_index(self):
        return pd.read_table(os.path.join(self.dir, 'profiles.index'), 
                             names=['model_number', 'priority', 'profile_number'], 
                             skiprows=1, sep='\s+')
    
    ### PROFILES
    def load_profile(self, profile_number):
        prof = pd.read_table(
            os.path.join(self.dir, 'profile' + str(profile_number) + '.data'), 
            skiprows=5, sep='\s+',
            usecols=self.usecols_profiles)
        return prof
    
    def get_profiles(self, parallel=True):
        if not parallel:
            return [self.load_profile(profile_number) 
                for profile_number in tqdm(self.index.profile_number)]
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            return list(tqdm(executor.map(self.load_profile, self.index.profile_number),
                                 total=len(self.index.profile_number)))
    
    ### GYRE FILES
    def load_gyre(self, profile_number):
        prof = gyre.load_gyre(
            os.path.join(self.dir, 'profile' + str(profile_number) + '.data.GYRE'))
        return prof
    
    def get_gyres(self, parallel=True):
        if not parallel:
            return [self.load_gyre(profile_number) 
                    for profile_number in tqdm(self.index.profile_number)]
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            return list(tqdm(executor.map(self.load_gyre, self.index.profile_number),
                                 total=len(self.index.profile_number)))
    
    ### FREQUENCIES
    def load_freq(self, profile_number):
        freq_path = os.path.join(self.dir, 'profile' + str(profile_number) + '-freqs.dat')
        freq = pd.read_table(freq_path, sep='\s+', skiprows=5) if os.path.isfile(freq_path) else None
        return freq

    def get_freqs(self, parallel=True):
        if not parallel:
            return [self.load_freq(profile_number) 
                    for profile_number in tqdm(self.index.profile_number)]
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            return list(tqdm(executor.map(self.load_freq, self.index.profile_number),
                                 total=len(self.index.profile_number)))


class Grid:
    def __init__(self, dir, load_history_extras=None, cpus=None, usecols_profiles=None, usecols_history=None):
        self.dir = dir
        self._df = None
        self.load_history_extras = load_history_extras
        self.cpus = os.cpu_count() if cpus is None else cpus
        self.usecols_profiles = usecols_profiles
        self.usecols_history = usecols_history
    
    def __repr__(self):
        return f"Grid with {len(self.df)} tracks\nColumns: " + ', '.join(list(self.df.columns)) + "\n"+\
            str(self.df.select_dtypes(include=np.number).apply(pd.Series.unique))
    
    @property
    def df(self):
        if self._df is None:
            self._df = self.load()
        return self._df
    
    def parse_filename(self, filename):
        parts = filename.split('-')
        parameters = {p.split('_')[0]: float(p.split('_')[1]) 
                      for p in parts}
        return parameters
    
    def process_directory(self, d, load_history_extras):
        if os.path.isdir(os.path.join(self.dir, d)):
            parameters = self.parse_filename(d)
            track = Track(os.path.join(self.dir, d), 
                          parameters=parameters, 
                          load_history_extras=load_history_extras,
                          usecols_profiles=self.usecols_profiles, 
                          usecols_history=self.usecols_history, 
                          cpus=self.cpus) 
            parameters['Track'] = track
            return parameters
    
    """
    def load(self):
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            data = list(tqdm(executor.map(lambda d: self.process_directory(d, self.load_history_extras), os.listdir(self.dir)),
                            total=len(os.listdir(self.dir))))
        data = [d for d in data if d is not None]  # Filter out any None results
        self._df = pd.DataFrame(data)
        for col in self._df.columns:
            if col not in ['Track']:
                self._df[col] = self._df[col].astype(float)
        return self._df
    """
    
    def load(self):
        self._df = self.process_directories(os.listdir(self.dir))
        return self._df
    
    def process_directories(self, dirs):
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            data = list(tqdm(executor.map(lambda d: self.process_directory(d, self.load_history_extras), dirs),
                            total=len(dirs)))
        data = [d for d in data if d is not None]  # Filter out any None results

        df = pd.DataFrame(data)
        for col in df.columns:
            if col not in ['Track']:
                df[col] = df[col].astype(float)

        return df

    def update(self):
        existing_dirs = set([os.path.dirname(track.dir) for track in self.df['Track']])
        all_dirs = set(os.listdir(self.dir))

        # Find new directories
        new_dirs = all_dirs - existing_dirs
        
        if new_dirs:
            new_df = self.process_directories(new_dirs)
            # Append to existing dataframe
            self._df = self.df.append(new_df, ignore_index=True)
        else:
            print("No new tracks found.")
    
    def filter(self, params):
        condition = True
        for param, values in params.items():
            if not isinstance(values, list):
                values = [values]
            condition = condition & self.df[param].isin(values)
        return self.df[condition]
