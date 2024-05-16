#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 04:44:24 2020

@author: pitagoras
"""

from joblib import load

prob_diff_values = {'0': 0.0, '1': 0.25, '2': 0.425, '3': 0.74, '4': 0.97}

class InteractionPredictor:
    def __init__(self, path, diff = 0.0):
        self.model = load(path)
        self.min_diff = diff
    
    def probs_to_bool(self, probs):
        return probs[:,0] + self.min_diff < probs[:,1]
    
    def predict(self, corr_values):
        assert len(corr_values) == 6
        probs = self.model.predict_proba([corr_values])
        classes = self.probs_to_bool(probs)
        return classes[0]
    
    def predict_many(self, corr_values_list, log = False):
        assert len(corr_values_list[0]) == 6
        probs = self.model.predict_proba(corr_values_list)
        if log:
            print("\n".join([str(x) for x in probs]))
        classes = self.probs_to_bool(probs)
        return classes
