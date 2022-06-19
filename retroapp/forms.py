from django import forms

from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _

from retroapp import constants


class TestForm(forms.Form):
    test_form_data = forms.CharField(help_text="Enter anything.")

    def clean_test_form_data(self):
        data = self.cleaned_data['test_form_data']

        print("This is the data. Nothing to validate.")
        print(data)

        # Remember to always return the cleaned data.
        return data


class SearchForm(forms.Form):
    smiles_string = forms.CharField(
        # help_text="Enter the SMILES string"
    )
    molecular_property = forms.ChoiceField(
        # help_text="Select property to query", 
        choices=enumerate(constants.MOLECULE_PROPERTIES.keys()),
    )
    molecular_property_min = forms.FloatField(
        # help_text="Enter the min value for the selected property"
    )
    molecular_property_max = forms.FloatField(
        # help_text="Enter the max value for the selected property"
    )
    sorting_mode = forms.ChoiceField(
        # help_text="Select the sorting method for the results",
        choices=enumerate(constants.SORTING_OPTIONS),
    )


    def clean_smiles_string(self):
        return self.cleaned_data['smiles_string']


    def clean_molecular_property(self):
        selected_property = self.fields['molecular_property'].choices[int(self.cleaned_data['molecular_property'])][1]

        if selected_property not in constants.MOLECULE_PROPERTIES.keys():
            raise ValidationError(_("Invalid property selected."))

        return selected_property
        
    def clean_molecular_property_min(self):
        mol_prop = self.cleaned_data['molecular_property']
        mol_prop_min = self.cleaned_data['molecular_property_min']

        if mol_prop_min < constants.MOLECULE_PROPERTIES[mol_prop]['min']:
            raise ValidationError(_(f"Invalid minimum value entered for the property. The value cannot be lower than {constants.MOLECULE_PROPERTIES[mol_prop]['min']}."))
        
        if mol_prop_min > constants.MOLECULE_PROPERTIES[mol_prop]['max']:
            raise ValidationError(_(f"Invalid minimum value entered for the property. The value cannot be greater than {constants.MOLECULE_PROPERTIES[mol_prop]['max']}."))
        
        return self.cleaned_data['molecular_property_min']

    def clean_molecular_property_max(self):
        mol_prop = self.cleaned_data['molecular_property']
        mol_prop_max = self.cleaned_data['molecular_property_max']

        if mol_prop_max > constants.MOLECULE_PROPERTIES[mol_prop]['max']:
            raise ValidationError(_(f"Invalid maximum value entered for the property. The value cannot be greater than {constants.MOLECULE_PROPERTIES[mol_prop]['max']}."))
        
        if mol_prop_max < constants.MOLECULE_PROPERTIES[mol_prop]['min']:
            raise ValidationError(_(f"Invalid maximum value entered for the property. The value cannot be lower than {constants.MOLECULE_PROPERTIES[mol_prop]['min']}."))
        
        return self.cleaned_data['molecular_property_max']

    def clean_sorting_mode(self):
        mode = self.fields['sorting_mode'].choices[int(self.cleaned_data['sorting_mode'])][1]
        if mode not in constants.SORTING_OPTIONS:
            raise ValidationError(_("Invalid sorting mode selected."))

        return mode