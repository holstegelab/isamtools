#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int cmd_index   (int argc, char **argv);
int cmd_import  (int argc, char **argv);
int cmd_check   (int argc, char **argv);
int cmd_sort    (int argc, char **argv);
int cmd_view    (int argc, char **argv);
int cmd_annotate(int argc, char **argv);
int cmd_update  (int argc, char **argv);
int cmd_rehint  (int argc, char **argv);
int cmd_audit_tags(int argc, char **argv);

static void usage(void) {
    fputs(
        "Usage: isamtools <command> [options]\n"
        "\n"
        "Commands:\n"
        "  index      Build store .isam (QNAME index) for a BAM/CRAM canonical store.\n"
        "  import     Import an alternate alignment into a .isam alignment view.\n"
        "  check      Verify BAM/CRAM/SAM content against a store or view.\n"
        "  sort       Build coordinate / name / interval indexes on a .isam.\n"
        "  view       Emit reconstructed SAM/BAM from a .isam alignment view.\n"
        "  annotate   Add read tags or named subsets to a .isam.\n"
        "  audit-tags Inspect observed aux tags and classify their scope.\n"
        "  update     Mutate metadata of a .isam.\n"
        "  rehint     Re-resolve alignment-view store_read_id values by QNAME.\n"
        "\n"
        "Run `isamtools <command> --help` for command-specific options.\n",
        stderr);
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(); return 2; }
    const char *cmd = argv[1];
    int sub_argc = argc - 1;
    char **sub_argv = argv + 1;

    if (!strcmp(cmd, "index"))    return cmd_index   (sub_argc, sub_argv);
    if (!strcmp(cmd, "import"))   return cmd_import  (sub_argc, sub_argv);
    if (!strcmp(cmd, "check"))    return cmd_check   (sub_argc, sub_argv);
    if (!strcmp(cmd, "sort"))     return cmd_sort    (sub_argc, sub_argv);
    if (!strcmp(cmd, "view"))     return cmd_view    (sub_argc, sub_argv);
    if (!strcmp(cmd, "annotate")) return cmd_annotate(sub_argc, sub_argv);
    if (!strcmp(cmd, "audit-tags")) return cmd_audit_tags(sub_argc, sub_argv);
    if (!strcmp(cmd, "update"))   return cmd_update  (sub_argc, sub_argv);
    if (!strcmp(cmd, "rehint"))   return cmd_rehint  (sub_argc, sub_argv);

    if (!strcmp(cmd, "-h") || !strcmp(cmd, "--help") || !strcmp(cmd, "help")) {
        usage(); return 0;
    }
    fprintf(stderr, "isamtools: unknown command '%s'\n", cmd);
    usage();
    return 2;
}
