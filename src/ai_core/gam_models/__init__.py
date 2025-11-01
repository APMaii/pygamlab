import importlib.resources as resources

def list_models():
    return [p.name for p in resources.files(__package__).iterdir() if p.suffix == ".gam_ai"]
