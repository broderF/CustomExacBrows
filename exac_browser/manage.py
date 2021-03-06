from flask.ext.script import Manager
from exac import app
import exac

manager = Manager(app)


@manager.command
def hello():
    print "hello"


@manager.command
def load_db():
    exac.load_db()


@manager.command
def load_base_coverage():
    exac.load_base_coverage()

@manager.command
def drop_variants():
    exac.drop_variants()

@manager.command
def load_exac_variants_file(filepath):
    exac.load_exac_variants_file(filepath)

@manager.command
def add_cohort(cohort_name):
    exac.addCohort(cohort_name)

@manager.command
def add_cohorts():
    exac.addCohorts()

@manager.command
def load_variants_file(filepath,cohort_name):
    exac.load_variants_file(filepath,cohort_name)

@manager.command
def update_annotation(filepath):
    exac.update_variant_annotation(filepath)

@manager.command
def load_analyse(cohort_name,analyse_name):
    analyse.load_analyse(cohort_name,analyse_name)

@manager.command
def dump_variants(filepath):
    exac.dump_all_variants(filepath)

@manager.command
def reload_variants():
    exac.load_variants_file()
    exac.load_mnps()
    exac.precalculate_metrics()

@manager.command
def load_meta_data(filepath):
    exac.loadMetaData(filepath)

@manager.command
def update(variant_file):
    exac.update(variant_file)

@manager.command
def load_gene_models():
    exac.load_gene_models()


@manager.command
def load_cnv_models():
    exac.load_cnv_models()


@manager.command
def load_cnv_genes():
    exac.load_cnv_genes()


@manager.command
def drop_cnv_genes():
    exac.drop_cnv_genes()


@manager.command
def load_dbsnp_file():
    exac.load_dbsnp_file()


@manager.command
def load_constraint_information():
    exac.load_constraint_information()


@manager.command
def load_mnps():
    exac.load_mnps()


@manager.command
def create_cache():
    exac.create_cache()


@manager.command
def precalculate_metrics():
    exac.precalculate_metrics()

if __name__ == "__main__":
    manager.run()

