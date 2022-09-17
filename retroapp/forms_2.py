from dataclasses import fields
from django import forms

from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _

import re

from retroapp.constants import MOLECULE_PROPERTIES, SORTING_OPTIONS, ASCENDING, DESCENDING, NO_SORT
from retroapp.models import QueryDB, QueryPropertyDB

from pprint import pprint


class FormQuery(forms.ModelForm):
    class Meta:
        model = QueryDB
        fields = ["Q_smiles", "Q_notes"]
        labels = {
            "Q_smiles": "SMILES String",
            "Q_notes": "Notes",
        }
        widgets = {
            "Q_notes": forms.Textarea(attrs={'rows':5, 'style': 'width: 100%'}),
        }

    def clean_Q_smiles(self):
        return self.cleaned_data['Q_smiles'].strip()

    def save(self, commit=True):
        # do something with self.cleaned_data['temp_id']
        return super(FormQuery, self).save(commit=commit)


class FormQueryProperty(forms.ModelForm):
    property_value_range = forms.CharField(
        label="Value Range",
        required=False,
        # help_text="Enter the min value for the selected property"
    )

    class Meta:
        model = QueryPropertyDB
        fields = [
            "Property_name", 
            "Sorting_mode",
        ]
        labels = {
            "Property_name": "Molecular Property",
            "Sorting_mode": "Sorting Mode",
        }

    def clean_Property_name(self):
        selected_property = self.fields['Property_name'].choices[int(self.cleaned_data['Property_name'])][0]    # choice index value

        return selected_property

        
    def clean_property_value_range(self):
        mol_prop = list(MOLECULE_PROPERTIES.keys())[self.cleaned_data['Property_name']]
        mol_prop_val_range = self.cleaned_data['property_value_range']
        mol_prop_val_range = mol_prop_val_range.replace(" ", "")

        if mol_prop_val_range:
            delim_search_res = re.search("[0-9]-", mol_prop_val_range)

            if delim_search_res is None:
                # CASE 1: Target
                target_val = float(mol_prop_val_range)

                if MOLECULE_PROPERTIES[mol_prop]['max'] < \
                    target_val < \
                    MOLECULE_PROPERTIES[mol_prop]['min']:
                    raise ValidationError(_(f"Invalid target value entered for the property. The value must be between {MOLECULE_PROPERTIES[mol_prop]['min']} and {MOLECULE_PROPERTIES[mol_prop]['max']}."))
            
            else:
                # CASE 2: Range
                min_val = float(mol_prop_val_range[:delim_search_res.span()[0]+1])
                max_val = float(mol_prop_val_range[delim_search_res.span()[1]:])

                if max_val < min_val:
                    raise ValidationError("Invalid Value Range. Range must be defined as start_val-end_val where start_val < end_val.")

                if MOLECULE_PROPERTIES[mol_prop]['max'] < \
                    min_val < \
                    MOLECULE_PROPERTIES[mol_prop]['min']:
                    raise ValidationError(_(f"Invalid minimum value entered for the property. The value must be between {MOLECULE_PROPERTIES[mol_prop]['min']} and {MOLECULE_PROPERTIES[mol_prop]['max']}."))

                if MOLECULE_PROPERTIES[mol_prop]['max'] < \
                    max_val < \
                    MOLECULE_PROPERTIES[mol_prop]['min']:
                    raise ValidationError(_(f"Invalid maximum value entered for the property. The value must be between {MOLECULE_PROPERTIES[mol_prop]['min']} and {MOLECULE_PROPERTIES[mol_prop]['max']}."))
                
                mol_prop_val_range = f"{min_val} - {max_val}"

        return mol_prop_val_range


    def clean_Sorting_mode(self):
        if (self.cleaned_data['Sorting_mode'] == "") or (self.cleaned_data['Sorting_mode'] is None):
            mode = NO_SORT
        else:
            mode = int(self.cleaned_data['Sorting_mode'])    # choice actual value
        
        print(f"Sorting {mode = }")

        return mode
    

    def save(self, commit=True):
        # do something with self.cleaned_data['temp_id']
        return super(FormQueryProperty, self).save(commit=commit)
