%module spammpack
%{
#include "spamm.h"
%}

char *
spamm_version ();

struct spamm_timer_t *
spamm_timer_new ();

void
spamm_timer_add_event (int event, struct spamm_timer_t *timer);

void
spamm_timer_delete (struct spamm_timer_t **timer);

void
spamm_timer_start (struct spamm_timer_t *timer);

void
spamm_timer_stop (struct spamm_timer_t *timer);

void
spamm_timer_info (const struct spamm_timer_t *timer, char *infostring, const int maxlength);

char *
spamm_timer_get_string (const struct spamm_timer_t *timer);
