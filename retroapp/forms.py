from django import forms

from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _

import re

from retroapp.constants import MOLECULE_PROPERTIES, SORTING_OPTIONS, ASCENDING, DESCENDING, NO_SORT


class SMILESForm(forms.Form):
    smiles_string = forms.CharField(
        label="SMILES String",
        # help_text="Enter the SMILES string"
    )
    notes = forms.CharField(
        label="Notes",
        widget=forms.Textarea(attrs={'rows':5, 'style': 'width: 100%'}),
        required=False,
    )

    def clean_smiles_string(self):
        return self.cleaned_data['smiles_string']


class PropertyForm(forms.Form):
    molecular_property = forms.ChoiceField(
        label="Molecular Property",
        # help_text="Select property to query", 
        choices=enumerate(MOLECULE_PROPERTIES.keys()),
        initial=list(MOLECULE_PROPERTIES.keys()).index('Cetane Number'),
    )
    property_value_range = forms.CharField(
        label="Value Range",
        required=False,
        # help_text="Enter the min value for the selected property"
    )
    sorting_mode = forms.ChoiceField(
        label="Sorting Mode",
        # help_text="Select the sorting method for the results",
        choices=enumerate(SORTING_OPTIONS),
        required=False,
        initial=SORTING_OPTIONS.index(DESCENDING),
    )

    def clean_molecular_property(self):
        selected_property = self.fields['molecular_property'].choices[int(self.cleaned_data['molecular_property'])][1]

        if selected_property not in MOLECULE_PROPERTIES.keys():
            raise ValidationError(_("Invalid property selected."))

        return selected_property

        
    def clean_property_value_range(self):
        mol_prop = self.cleaned_data['molecular_property']
        # mol_prop_min = self.cleaned_data['molecular_property_min']
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


    def clean_sorting_mode(self):
        if self.cleaned_data['sorting_mode'] == "":
            mode = NO_SORT
        else:
            mode = self.fields['sorting_mode'].choices[int(self.cleaned_data['sorting_mode'])][1]
        
        if mode not in SORTING_OPTIONS:
            raise ValidationError(_("Invalid sorting mode selected."))

        return mode