from mp_api.client import MPRester

# Option 1: Pass your API key directly as an argument.
with MPRester("YOUR api") as mpr:
    # user stuff with mpr...
    print(mpr.get_structure_by_material_id("mp-4019").formula)
# Option 2: Use the `MP_API_KEY` environment variable:
# export MP_API_KEY="your_api_key_here"
# Note: You can also configure your API key through pymatgen